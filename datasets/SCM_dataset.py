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

from .data_wrapper import DataWrapper


class SCM_dataset(DataWrapper):
    def __init__(self, cfg, logger=None):
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

        self.data = np.load(self.data_path, allow_pickle=True)
        flag= False
        for item in self.data:
            if item['idx'] == self.idx:
                self.data = item
                flag = True
                break
        if not flag:
            raise ValueError(f"Index {self.idx} not found in the dataset at {self.data_path}.")

        self.time_series = self.data['data']
        self.true_edges = self.data['true_edges']
        self.graph_str = self.draw_graph(self.true_edges)
        self.num_vars = self.time_series.shape[-1]
        self.dim = self.num_vars
        self.var_list = [f'x_{i}' for i in range(self.num_vars)]

        if 'synthetic' in self.data_path.lower():
            self.process_syntetic()
        elif 'science'in self.data_path.lower():
            self.process_science()
        else:
            raise ValueError("Unknown dataset type. Please specify 'strogatz' or 'physiome' in the data path.")
        
        if self.cfg.experiment.symbolic_expression.get("use_context", False) and self.cfg.experiment.symbolic_expression.context.get("image", False): 
            # self.image = f"data/Synthetic_SCM_images/{self.idx}.png"
            pass
        else:
            self.image = None
    

    def draw_graph(self, edge_index):
        """
        Convert edge_index list of (src, dst, lag) tuples into a readable edge list string.
        """
        # each tuple is an individual edge with a lag value
        if type(edge_index) is list:
            lines = [f"{src} -> {dst} (lag={lag})" for src, dst, lag in edge_index]
            graph_str = " | ".join(lines)
            if self.logger:
                self.logger.info("Graph edges with lag:\n" + graph_str)
            else:
                print(graph_str)
            return graph_str
        elif isinstance(edge_index, dict):
            # group edges by target: 'target <- src1_(lag), src2_(lag), ...'
            lines = []
            for dst in sorted(edge_index):
                parents = edge_index[dst]
                # format each parent as src_(lag=...)
                parts = [f"{src}_(lag={lag})" for src, lag in parents]
                lines.append(f"{dst} <- " + ", ".join(parts))
            graph_str = " | ".join(lines)
            if self.logger:
                self.logger.info("Graph edges with lag:\n" + graph_str)
            else:
                print(graph_str)
            return graph_str
         
    def process_synthetic(self):
         # import ipdb; ipdb.set_trace()
        # self.eq_description = self.data['eq_description'].item()
        self.lagging = self.data['L_max']

        self.var_description = None
        self.domain = None
        
        self.logger.info(f"Loaded data from {self.data_path} with index {self.idx}.")
        self.logger.info(f"Groundtruth: {self.data['true_edges']}")
    
    def process_science(self):
        # import ipdb; ipdb.set_trace()
        self.eqs_str = self.data['eq_str']
        self.var_description = str(self.data['variable_names'])
        self.domain = self.data['domain']
        self.lagging = 1
        
        self.logger.info(f"Loaded data from {self.data_path} with index {self.idx}.")
        self.logger.info(f"Groundtruth: {self.data['true_edges']}")


