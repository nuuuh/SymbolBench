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
from sklearn.metrics import precision_score, recall_score, f1_score

class DataWrapper:
    def __init__(self, cfg):
        self.cfg = cfg
    
    

    


