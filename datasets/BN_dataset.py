import torch
import numpy as np
import pandas as pd
import sympy
import ipdb
from tqdm import tqdm
import json
import os

from odeformer.odebench.strogatz_equations import equations
from odeformer.odebench.solve_and_plot import config, process_equations, solve_equations, plot_prediction

from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from sympy import count_ops
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score

from .data_wrapper import DataWrapper



class BN_dataset(DataWrapper):
    def __init__(self, cfg, logger=None):
        super().__init__(cfg)
        self.cfg = cfg
        # self.logger = logger
        self.root_dir = cfg.get("root", os.getcwd())

        self.idx = cfg.experiment.symbolic_expression.get("idx", None)
        if self.idx is None:
            raise ValueError("Index is not specified in the configuration.")
        self.data_path = cfg.experiment.symbolic_expression.get("data_path", None)
        self.data_path = os.path.join(self.root_dir, self.data_path) if self.data_path is not None else None
        if self.data_path is None:
            raise ValueError("Data path is not specified in the configuration.")
        try:
            with open(self.data_path, 'r') as f:
                self.data = json.load(f)
            self.data = pd.DataFrame(self.data)
            self.data = self.data[self.data['idx'] == self.idx].iloc[0]
        except:
            try:
                print(f"Loading data from {self.data_path} failed. Attempting to load from CSV.")
                self.data = pd.read_csv("data/boolnet.csv")
                # import ipdb; ipdb.set_trace()
                input_obs_list = []
                transitions = []
                for j in range(30):
                    bool_data = self.get_data(self.idx)
                    data = bool_data['sim'].T.to_dict()
                    input_obs = bool_data['input_observations']
                    timepoints = sorted(t for t in data if isinstance(t, int))
                    trans = [(data[t], data[t+1]) for t in timepoints[:-1]]
                    transitions += [ t for t in trans if t not in transitions ]
                    input_obs_list.append(input_obs)
                
                train_transitions = transitions[:len(transitions)//2]
                test_transitions = transitions[len(transitions)//2:]
                
                self.data = self.data[self.data['idx'] == self.idx].iloc[0]

                self.data['train_transitions'] = train_transitions
                self.data['test_transitions'] = test_transitions
                # import ipdb; ipdb.set_trace()
                self.data['initial_obs'] = [obs.to_dict() for obs in input_obs_list]

                # Save updated data: load existing JSON, replace or append the entry matching self.idx
                if os.path.exists(self.data_path):
                    try:
                        with open(self.data_path, 'r') as fj:
                            data_list = json.load(fj)
                    except (json.JSONDecodeError, ValueError):
                        # file empty or invalid JSON, start fresh list
                        data_list = []
                else:
                    data_list = []
                
                updated_entry = self.data.to_dict() if hasattr(self.data, 'to_dict') else self.data
                for i, item in enumerate(data_list):
                    if item.get('id') == self.idx:
                        data_list[i] = updated_entry
                        break
                else:
                    data_list.append(updated_entry)

                with open(self.data_path, 'w') as f:
                    json.dump(data_list, f, indent=2)

            except Exception as e:
                raise ValueError(f"Failed to load data from {self.data_path}. Error: {e}")


        self.transitions_str = self.transitions_to_string(self.data['train_transitions'])

        self.context = self.data['model_name']
        self.num_vars = self.data['num_var']
        self.num_eqs = self.data['num_regulars']
        self.variable_mapping = self.data['variable_mapping']
        self.gt_exprs_str = self.data['regulations']

        self.trans_vars, self.gt_exprs, _  = self.parse_regulations_to_sympy(self.gt_exprs_str)
        self.var_list = [v.split(':')[0].strip() for v in self.variable_mapping.split(',')]
        self.free_vars = [v for v in self.var_list if v not in self.trans_vars]

        # import ipdb; ipdb.set_trace()

    def transitions_to_string(self, transition_list):
        str_transitions = []
        for transitions in transition_list:
            before = transitions[0]
            after = transitions[1]

            variables = sorted(before.keys(), key=lambda x: int(x[1:]))  # ['x1','x2',…]
            before_str = ''.join('1' if before[v] else '0' for v in variables)
            after_str  = ''.join('1' if after[v]  else '0' for v in variables)
            str_transitions.append(f"{before_str} -> {after_str}")
        return '; '.join(str_transitions)

    def get_observed_inputs(self, sim_length, variables, inputs, init_state=None, observations=None):
        """
        Generate observed inputs for the simulation.
        """
        if init_state is None:
            init_state = {v: bool(np.random.randint(2)) for v in variables}
        if observations is None:
            observations = pd.DataFrame(
                [[bool(np.random.randint(2)) for j in range(len(variables))] for i in range(sim_length)],
                columns=variables
            )

        rows = [init_state]
        for i in range(1, sim_length):
            row_dict = {inp: observations.iloc[i][inp] for inp in inputs}
            rows.append(row_dict)
        input_observations = pd.DataFrame(rows, columns=variables)

        return observations, input_observations
    
    def get_all_data(self, sim_length=10):
        indices = self.data.index
        results = []
        
        for idx in tqdm(indices):
            row = self.data.iloc[idx]
            if sim_length is None:
                sim_length = row['simulation_length']

            variables = [v.strip() for v in row['variable_mapping'].split(',')]
            variables = [v.split(':')[0] for v in variables]
            reg_lines = row['regulations'].split('\n')
            reg_vars = [line.split('=')[0].strip() for line in reg_lines if '=' in line]
            inputs = list(tuple(set(variables) - set(reg_vars)))

            input_obs, input_observations = self.get_observed_inputs(sim_length, variables, inputs, init_state=None, observations=None)
            # import ipdb; ipdb.set_trace()

            row = self.data.iloc[idx]

            truth_tables, expr, attractor_info = self.generate_boolean_network_data(
                                        simulation_length=sim_length,
                                        regulations = row['regulations'],
                                        input_observations=input_observations,
                                        mapping=row['variable_mapping'],
                                    )

            row_result = {"sim": truth_tables,
                            "inputs": inputs, 
                            "vars": reg_vars, 
                            "model_index": idx, 
                            "input_observations": input_observations, 
                            "regulations": expr, 
                            "regulations_str": row['regulations'],
                            "observations": input_obs, 
                            "mapping": row['variable_mapping'], 
                            "attractor_info": attractor_info}                     
            results.append(row_result)

        return results

    def get_data(self, idx, sim_length=10):
        """
        Get the data for a specific index.
        """
        if idx is None:
            return self.get_all_data(sim_length=sim_length)
        else:
            row = self.data[self.data.idx == idx].iloc[0]
            if sim_length is None:
                sim_length = row['simulation_length']

            variables = [v.strip() for v in row['variable_mapping'].split(',')]
            variables = [v.split(':')[0] for v in variables]
            reg_lines = row['regulations'].split('\n')
            reg_vars = [line.split('=')[0].strip() for line in reg_lines if '=' in line]
            inputs = list(tuple(set(variables) - set(reg_vars)))

            input_obs, input_observations = self.get_observed_inputs(sim_length, variables, inputs, init_state=None, observations=None)

            truth_tables, expr, attractor_info = self.generate_boolean_network_data(
                                        simulation_length=sim_length,
                                        regulations = row['regulations'],
                                        input_observations=input_observations,
                                        mapping=row['variable_mapping'],
                                    )

            row_result = {"sim": truth_tables,
                            "inputs": inputs, 
                            "vars": reg_vars, 
                            "model_index": idx, 
                            "input_observations": input_observations, 
                            "regulations": expr, 
                            "regulations_str": row['regulations'],
                            "observations": input_obs, 
                            "mapping": row['variable_mapping'], 
                            "attractor_info": attractor_info}                     

            return row_result

            

    def detect_attractor(self, df):
        """
        Detects if the simulation reaches a stable state (attractor).
        Returns (is_attractor_found, attractor_start, attractor_length, attractor_states)
        """
        seen = {}
        for t, row in df.iterrows():
            state_tuple = tuple(row.values)
            if state_tuple in seen:
                attractor_start = seen[state_tuple]
                attractor_length = t - attractor_start
                attractor_states = df.iloc[attractor_start:t]
                return True, attractor_start, attractor_length, attractor_states
            seen[state_tuple] = t
        return False, None, None, None

    def generate_boolean_network_data(
        self,
        regulations=None,
        simulation_length=None,
        input_observations=None,
        mapping=None,
    ):
        """
        Generate boolean network data for a single model (specified by model_index).
        """
        # import ipdb; ipdb.set_trace()
        if type(regulations) == str:
            variables, exprs, free_vars = self.parse_regulations_to_sympy(regulations)
        elif type(regulations) == dict:
            variables = list(tuple(regulations.keys()))
            # import ipdb; ipdb.set_trace()
            exprs = regulations
        else:
            raise ValueError("Regulations must be a string or a dictionary.")

        if mapping is not None:
            all_vars = [v.split(':')[0].strip() for v in mapping.split(',')]
        else:
            all_vars = input_observations.columns.to_list()

        # reg_vars = [line.split('=')[0].strip() for line in regulations_str.split("\n") if '=' in line]
        # ipdb.set_trace()
        input_vars = list(tuple(set(all_vars) - set(variables)))


        # --- Simulation mode ---
        sim_obs = input_observations[input_vars].copy()
        init_state = input_observations.iloc[0][all_vars].to_dict()

        states = [init_state.copy()]

        # ipdb.set_trace()
        for t in range(1, simulation_length):

            prev_state = states[-1].copy()
            # Ensure all variables have bool values (no NaN)
            for v in all_vars:
                if v not in prev_state or pd.isna(prev_state[v]):
                    prev_state[v] = False
            # Compute new state synchronously from previous state
            # ipdb.set_trace()
            next_state = prev_state.copy()
            for v, expr in exprs.items():
            # for v, expr in zip(variables, exprs):
                next_state[v] = bool(expr.subs(prev_state))
            states.append(next_state)

        df = pd.DataFrame(states, columns=all_vars)
        df.index.name = 'time_step'

        # Attractor detection
        found, start, length, attractor = self.detect_attractor(df)
        attractor_info = {
            'found': found,
            'start': start,
            'length': length,
            'states': attractor
        }
        # Truncate output to only include until attractor (including first attractor state)
        if found:
            df = df.iloc[:start+length]
        return df, exprs, attractor_info

    @staticmethod
    def parse_regulations_to_sympy(regulations_str):
        exprs = {}
        variables = []
        free_vars = []
        from sympy.parsing.sympy_parser import parse_expr, standard_transformations
        transformations = standard_transformations
        # import ipdb; ipdb.set_trace()
        if ';' not in regulations_str:
            eqs = regulations_str.strip().split('\n')
        else:
            eqs = regulations_str.strip().split(';')
        for line in eqs:
            if '=' in line:
                lhs, rhs = line.split('=', 1)
                lhs = lhs.strip()
                rhs = rhs.strip()
                rhs = rhs.replace('AND', '&').replace('XOR', '^').replace('OR', '|').replace('NOT', '~').replace('!', '~')
                expr = parse_expr(rhs, transformations=transformations)
                exprs[lhs] = expr
                variables.append(lhs)
                free_vars.append(sorted(str(v) for v in expr.free_symbols))
            
        return variables, exprs, free_vars

    @staticmethod
    def get_global_free_variables(regulations_str):
        variables, exprs, free_vars = BN_dataset.parse_regulations_to_sympy(regulations_str)
        all_free_vars = set(var for sublist in free_vars for var in sublist)
        dependent_vars = set(variables)
        global_free_vars = sorted(all_free_vars - dependent_vars)
        return global_free_vars
    