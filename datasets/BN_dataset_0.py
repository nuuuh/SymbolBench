import torch
import numpy as np
import pandas as pd
import sympy
import ipdb
from tqdm import tqdm
import json

from odeformer.odebench.strogatz_equations import equations
from odeformer.odebench.solve_and_plot import config, process_equations, solve_equations, plot_prediction

from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from sympy import count_ops
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score

from .data_wrapper import DataWrapper


class BN_dataset_0(DataWrapper):
    def __init__(self, data):
        super().__init__(data)
        self.data = data

    def print_data(self, n=5):
        print(self.data.head(n))

    def __str__(self):
        return f"BN_dataset with {len(self.data)} entries. Columns: {list(self.data.columns)}"

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
            row_result = {
                "idx": row.idx,
                "sim": truth_tables,
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
            row = self.data.iloc[idx]
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
        for line in regulations_str.strip().split('\n'):
            if '=' in line:
                lhs, rhs = line.split('=', 1)
                lhs = lhs.strip()
                rhs = rhs.strip()
                rhs = rhs.replace('AND', '&').replace('OR', '|').replace('NOT', '~').replace('!', '~')
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

    def evaluate_model(self, gt_exprs, inferred, vairables, inputs, input_observations, mapping, ood=False):
        variables = sorted(vairables, key=lambda x: int(x.split('x')[1]))
        # evaluate complexity

        # evaluate accuracy
        # evaluate precision, recall, f1
        
        if ood:
            new_input_observations = input_observations.copy()
            while new_input_observations.equals(input_observations):
                new_input_obs, new_input_observations = self.get_observed_inputs(len(input_observations), vairables, inputs, init_state=None, observations=None)
            input_observations = new_input_observations
        
        gt_data = self.generate_boolean_network_data(
                                                simulation_length=len(input_observations),
                                                regulations = gt_exprs,
                                                input_observations=input_observations,
                                                mapping=mapping)
        
        gt_solutions = gt_data[0]
        num_nodes = len(gt_exprs)
        gt_ops = sum(expr.count_ops() for expr in gt_exprs.values())
        avg_gt_ops = gt_ops / num_nodes

        try:
            pred_data = self.generate_boolean_network_data(
                                                    simulation_length=len(input_observations),
                                                    regulations = inferred,
                                                    input_observations=input_observations,
                                                    mapping=mapping)
            pred_solutions = pred_data[0]
            pred_ops = sum(expr.count_ops() for expr in inferred.values())
            avg_pred_ops = pred_ops / num_nodes
        except Exception as e:
            print(f"Error in generating predicted data: {e}")

            return {
                    'gt_data': gt_solutions,
                    'pred_data': None,
                    'gt_complexity': avg_gt_ops,
                    'pred_complexity': avg_gt_ops,
                    'precision': 0,
                    'recall': 0,
                    'f1': 0,
                    'accuracy': 0,
                    }
        
        
        # import ipdb; ipdb.set_trace()
        while len(gt_solutions) > len(pred_solutions):
            pred_solutions = pred_solutions._append(pred_solutions.iloc[-1])

        while len(gt_solutions) < len(pred_solutions):
            gt_solutions = gt_solutions._append(gt_solutions.iloc[-1])
        
        pred_solutions.index = range(len(pred_solutions))
        gt_solutions.index = range(len(gt_solutions))

        # compute precision, recall, f1
        y_true = gt_solutions.values.flatten()
        y_pred = pred_solutions.values.flatten()

        precision = precision_score(y_true, y_pred)
        recall = recall_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)
        # compute accuracy
        accuracy = accuracy_score(y_true, y_pred)

        return {
            'gt_data': gt_solutions,
            'pred_data': pred_solutions,
            'gt_complexity': avg_gt_ops,
            'pred_complexity': avg_pred_ops,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'accuracy': accuracy,
        }