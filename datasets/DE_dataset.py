import torch
import numpy as np
import pandas as pd
import sympy
import ipdb
from tqdm import tqdm
import json
import os
import re

from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from sympy import count_ops
from sklearn.metrics import precision_score, recall_score, f1_score

from odeformer.odebench.strogatz_equations import equations
from odeformer.odebench.solve_and_plot import config, process_equations, solve_equations, plot_prediction
from .data_wrapper import DataWrapper


class DE_dataset(DataWrapper):
    def __init__(self, cfg, logger):
        super().__init__(cfg)
        self.cfg = cfg
        self.logger = logger
        self.root_dir = cfg.get("root", os.getcwd())
        self.idx = cfg.experiment.symbolic_expression.get("idx", None)
        if self.idx is None:
            raise ValueError("Index is not specified in the configuration.")
        self.data_path = cfg.experiment.symbolic_expression.get("data_path", None)
        self.data_path = os.path.join(self.root_dir, self.data_path) if self.data_path is not None else None
        if self.data_path is None:
            raise ValueError("Data path is not specified in the configuration.")
        # import ipdb; ipdb.set_trace()
        with open(self.data_path, 'r') as f:
            self.data = json.load(f)
        self.data = pd.DataFrame(self.data)
        
        self.data = self.data[self.data['id'] == self.idx]

        if 'strogatz' in self.data_path.lower():
            self.process_strogatz()
        elif 'physiome' in self.data_path.lower():
            self.process_physiome()
        elif "DE_2_dims" in self.data_path.lower():
            self.process_strogatz()
        else:
            self.process_analysis()
            logger.info("Unknown dataset type. Using analysis dataset.")
        


    def __get_num_vars__(self):
        pass

    def __get_num_eqs__(self):
        pass
    
    def get_data(self, idx):
        return self.data.iloc[idx]
        
    def str_to_sympy(self, eq_string, var_list=None, consts=None):
        allowed_operators = set(['+', '-', '*', '/', '**', '^'])
        allowed_functions = set([
            'sqrt', 'exp', 'log', 'abs',
            'sin', 'cos', 'cot', 'tan', 'sinh', 'cosh', 'tanh'
        ])

        individual_eqs = eq_string.split('|')
        parsed_eqs = []
        for eq in individual_eqs:
            # import ipdb; ipdb.set_trace()
            expr = sympy.sympify(eq, evaluate=True)
            # Check for allowed functions (case-insensitive)
            for func in expr.atoms(sympy.Function):
                fname = func.func.__name__.lower()
                if fname not in allowed_functions:
                    raise ValueError(f"Disallowed function: {func.func.__name__} in {eq}")
            # Check for allowed operators
            for node in sympy.preorder_traversal(expr):
                if isinstance(node, sympy.Pow):
                    op = '^'
                elif isinstance(node, sympy.Mul):
                    op = '*'
                elif isinstance(node, sympy.Add):
                    op = '+'
                # No need to check for Sub or Div
                else:
                    continue
                if op not in allowed_operators:
                    raise ValueError(f"Disallowed operator: {op} in {eq}")
                
            parsed_eqs.append(expr)
        
        if consts is not None:
            # import ipdb; ipdb.set_trace()
            extracted_consts = re.findall(r'c_\d+', eq_string)
            consts_values = [consts[int(const.split('_')[-1])] for const in extracted_consts]
            const_symbols = sympy.symbols(extracted_consts)
            const_subs = dict(zip(const_symbols, consts_values))
            # import ipdb; ipdb.set_trace()

            return [eq.subs(const_subs) for eq in parsed_eqs]
        else:
            return parsed_eqs

    def process_strogatz(self):
        # import ipdb; ipdb.set_trace()
        self.eqs_str = self.data['eq'].item()
        self.consts = self.data['consts'].item()
        self.initial_conditions = self.data.init.item()[0] # the first initial condition
        # self.eq_description = self.data['eq_description'].item()
        self.const_desciption = self.data['const_description'].item()
        self.var_description = self.data['var_description'].item()
        self.consts_range = self.data.const_constraints.item()

        # import ipdb; ipdb.set_trace()
        self.num_eqs = int(self.data.dim.item())
        self.eqs = self.str_to_sympy(self.eqs_str, None, self.consts[0])
        self.num_vars = max(len(set().union(*[eq.free_symbols for eq in self.eqs])), self.num_eqs)
        self.var_symbols = sympy.symbols([f'x_{i}' for i in range(self.num_vars)])
       
        self.image = f"data/Strogatz_images/{self.idx}.png"
        self.domain = None
        
        self.logger.info(f"Loaded data from {self.data_path} with index {self.idx}.")
        
        self.logger.info(f"Groundtruth Function: {self.data}")
    
    def process_physiome(self):
        # import ipdb; ipdb.set_trace()
        self.eqs_str = self.data['equation'].item()
        self.consts = self.data['constants'].item()
        self.initial_conditions = self.data.initial_states.item()
        self.const_desciption = str(self.data['c_mapping'].item())
        self.var_description = str(self.data['x_mapping'].item())
        self.consts_range = None

        # import ipdb; ipdb.set_trace()
        self.num_eqs = int(self.data.dim.item())
        self.eqs = self.str_to_sympy(self.eqs_str, None, self.consts)
        self.num_vars = len(set().union(*[eq.free_symbols for eq in self.eqs]))
        self.var_symbols = sympy.symbols([f'x_{i}' for i in range(self.num_eqs)])

        self.image = f"data/Physiome_images/{self.idx}.png"
        self.domain = self.data['domain'].item()
        
        self.logger.info(f"Loaded data from {self.data_path} with index {self.idx}.")
        
        self.logger.info(f"Groundtruth Function: {self.data}")

    def process_analysis(self):
        self.eqs_str = self.data['equation'].item()
        self.consts = self.data['constants'].item()
        self.initial_conditions = self.data.initial_states.item()
        self.const_desciption = str(self.data['c_mapping'].item())
        self.var_description = str(self.data['x_mapping'].item())
        self.consts_range = None

        # import ipdb; ipdb.set_trace()
        self.num_eqs = int(self.data.dim.item())
        self.eqs = self.str_to_sympy(self.eqs_str, None, self.consts)
        self.num_vars = len(set().union(*[eq.free_symbols for eq in self.eqs]))
        self.var_symbols = sympy.symbols([f'x_{i}' for i in range(self.num_eqs)])

        self.image = None
        try:
            self.domain = self.data['domain'].item()
        except:
            self.domain = None
        
        self.logger.info(f"Loaded data from {self.data_path} with index {self.idx}.")
        
        self.logger.info(f"Groundtruth Function: {self.data}")
    
        
        
