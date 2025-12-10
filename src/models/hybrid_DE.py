import os
from utils import utils
import numpy as np
import random
from utils.data_utils import ExpressionTree

from gplearn.genetic import SymbolicRegressor
from gplearn.fitness import make_fitness
import sklearn
from sklearn.utils.validation import validate_data
if not hasattr(SymbolicRegressor, '_validate_data'):
    SymbolicRegressor._validate_data = validate_data

import torch
import re
from odeformer.envs.generators import NodeList
import sympy
import sympy as sp
from scipy.integrate import solve_ivp
from sympy import sympify, lambdify, symbols, sin, cos, log
from transformers import AutoTokenizer, AutoModelForCausalLM
# Assuming TSProcessor and snapshot_download are available in the environment
# from some_module import TSProcessor, snapshot_download
from odeformer.utils import timeout, MyTimeoutError
import pickle
import tempfile
from odeformer.metrics import r2_score
from utils.utils import compute_ned

from .hf_model import HuggingFaceModel
from .openai_model import OpenAIModel


class Hybrid_DE(object):
    def __init__(self, model_name, device, dtype, cache_dir=None, **kwargs):
        self.set_seed(15)
        self.model_name = model_name
        self.device = device
        self.dtype = dtype
        token = os.environ.get("HF_TOKEN", None)

        if "gpt" not in model_name:
            self.llm = HuggingFaceModel(model_name, device, dtype, cache_dir, **kwargs)
        elif "gpt" in model_name:
            self.llm = OpenAIModel(model_name, device, dtype, cache_dir, **kwargs)
        self.simplicity_weight = kwargs.get("simplicity_weight", 0.05)
        
        # Parsimony coefficient for gplearn
        self.complexity_weight = kwargs.get("complexity_weight", 0.1)

        # Initialize gplearn SymbolicRegressor for each variable
        gp_kwargs = {
            "population_size": 25,
            "generations": 100,
            "function_set": ['add', 'sub', 'mul', 'div', 'sin', 'cos', 'log'],
            "parsimony_coefficient": self.complexity_weight,
            "p_crossover": 0.9,
            "p_subtree_mutation": 0.01,
            "p_hoist_mutation": 0.01,
            "p_point_mutation": 0.01,
            "max_samples": 1.0,
            "verbose": 1,
            "random_state": 0,
        }
        self.gp_kwargs = gp_kwargs
        self.gp_models = []  # one SymbolicRegressor per variable

        self.use_llm=False
        self.llm_weight=1
        # import ipdb; ipdb.set_trace()
        custom_metric = make_fitness(function=self.symbolic_mse, greater_is_better=False)
        # include custom metric into gp parameters
        self.gp_kwargs['metric'] = custom_metric

        self.solve_config = {
            "t_span": (0, 10),
            "method": "LSODA",       # switch to BDF to avoid LSODA prints
            "rtol": 1e-6,
            "atol": 1e-6,
            "first_step": 1e-6,
            "t_eval": np.linspace(0, 10, 150),
            "dense_output": False,
        }
    

    def set_seed(self, seed):
        """Set random seed for reproducibility."""
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)


    def gplearn_to_sympy(self, expr_str):
        # collect all X-variable names (X0, X1, …)
        try:
            var_names = sorted(set(re.findall(r'X\d+', expr_str)),
                            key=lambda s: int(s[1:]))

            # otherwise build symbols and mapping
            # syms = symbols(' '.join(var_names))
            syms = symbols(var_names)
            local_dict = dict(zip(var_names, syms))

            # map gplearn fun names to Sympy operators
            local_dict.update({
                'add': lambda x, y: x + y,
                'sub': lambda x, y: x - y,
                'mul': lambda x, y: x * y,
                'div': lambda x, y: x / y,
                'sin': sin,
                'cos': cos,
                'log': log
            })

            simplified = sympify(expr_str, locals=local_dict)
        except Exception as e:
            print(f"Error in gplearn_to_sympy: {e}")
            simplified = expr_str


        return simplified


    def symbolic_mse(self, y, y_pred, weights, eq_str):
        # import ipdb; ipdb.set_trace()
        eq_str = self.gplearn_to_sympy(eq_str)

        mse_loss = np.mean((y - y_pred)**2)

        llm_score=0

        if self.use_llm:
            # import ipdb; ipdb.set_trace()
            scoring_prompt = self.prompt.replace("{candidate_exprs}", str(eq_str))
            llm_response = self.llm.generate(scoring_prompt)
            try:
                llm_score = float(re.findall(r"[-+]?\d*\.\d+|\d+", llm_response)[0])
            except (ValueError, IndexError):
                llm_score = 0.0

            llm_score = self.llm_weight*(1- llm_score/10)
        if type(eq_str) is float or type(eq_str) is int:
            complexity = 0
        else:
            complexity = eq_str.count_ops()

        return mse_loss + llm_score + 0*complexity


    def fit(self, trajectories, times=None, **kwargs):
        # Compute time points
        self.current_trajectories = trajectories
        # import ipdb; ipdb.set_trace()
        # Estimate derivatives via finite differences
        dt = np.gradient(times)
        derivatives = np.gradient(trajectories, axis=0) / dt[:, None]
        # Prepare training inputs
        X = trajectories
        # Fit one SymbolicRegressor per variable using stored gp_kwargs
        
        for i in range(trajectories.shape[1]):
            model = SymbolicRegressor(**self.gp_kwargs)
            model.fit(X, derivatives[:, i])
            self.gp_models.append(model)
        

        fitted_exprs = [re.sub(r'X(\d+)', r'x_\1', str(item)) for item in self.gp_models]
        fitted_exprs = [self.gplearn_to_sympy(expr) for expr in fitted_exprs]

        # import ipdb; ipdb.set_trace()
        return fitted_exprs

    def get_best_equations(self):
        # Return expressions for each variable
        equations = {i: [model._program.__str__()] for i, model in enumerate(self.gp_models)}
        return equations

    def predict(self, times, initial_conditions=None):
        """Make predictions using the best found equations."""
        # Build ODE system from learned expressions
        # Prepare lambdified functions
        funcs = []
        for model in self.gp_models:
            expr = sympify(model._program.__str__())
            f = lambdify([*['x'+str(i) for i in range(len(self.gp_models))]], expr, 'numpy')
            funcs.append(f)
        # Define ODE rhs
        def rhs(t, y):
            return [func(*y) for func in funcs]
        # Set initial conditions
        y0 = initial_conditions if initial_conditions is not None else self.current_trajectories[0]
        sol = solve_ivp(rhs, (times[0], times[-1]), y0, t_eval=times)
        return sol.y.T
    

    def final_evaluate(self, inferred, dataset):
        idx = dataset.idx
        gt_func = dataset.eqs
        try:
            gt_exprs = gt_func
            inferred_exprs = inferred

            dim = len(gt_exprs)
            var_symbols = sp.symbols([f'x_{i}' for i in range(dim)])

            # Ensure constants (e.g., pi, E) are numeric floats
            gt_exprs = [expr.evalf() for expr in gt_exprs]
            inferred_exprs = [expr.evalf() for expr in inferred_exprs]

            gt_funcs = [sp.lambdify(var_symbols, expr, 'numpy') for expr in gt_exprs]
            inferred_funcs = [sp.lambdify(var_symbols, expr, 'numpy') for expr in inferred_exprs]

            # Define ODE functions
            def gt_ode(t, x):
                return np.array([f(*x) for f in gt_funcs])

            def inferred_ode(t, x):
                return np.array([f(*x) for f in inferred_funcs])

            # Solve ODEs
            # ID
            if type(dataset.data.init.item()[0]) is list and len(dataset.data.init.item()) > 0:
                ID_initial_conditions = dataset.data.init.item()[0]
                OOD_initial_conditions = dataset.data.init.item()[1]
            else:
                ID_initial_conditions = dataset.initial_conditions
                OOD_initial_conditions = np.abs(ID_initial_conditions + np.random.normal(0, 0.1, size=len(ID_initial_conditions)))

            def generate(init_conds):
                gt_sol = solve_ivp(gt_ode, **self.solve_config, y0=init_conds)
                inferred_sol = solve_ivp(inferred_ode, **self.solve_config, y0=init_conds)
                gt_trajectories = gt_sol.y.T
                inferred_trajectories = inferred_sol.y.T
                R2 = r2_score(gt_trajectories, inferred_trajectories)
                return R2

            # import ipdb; ipdb.set_trace()
            ID_R2 = generate(ID_initial_conditions)
            OOD_R2 = generate(OOD_initial_conditions)

            gt_complexity = sum(expr.count_ops() for expr in gt_exprs)
            inferred_complexity = sum(expr.count_ops() for expr in inferred_exprs)
            ned_metrics = sum(compute_ned(inferred_exprs, gt_exprs, var_symbols))
            results = {
                'index': idx,
                'num_eqs': dim,
                'predicted_equation': inferred_exprs,
                'groundtruth_equation': gt_exprs,
                'ID_R2': ID_R2,
                'OOD_R2': OOD_R2,
                'predicted_complexity': inferred_complexity,
                'groundtruth_complexity': gt_complexity,
                'initial_conditions': [ID_initial_conditions, OOD_initial_conditions],
                'ned': ned_metrics,
                }
        

        except Exception as e:
            results = {
                'index': idx,
                'num_eqs': None,
                'predicted_equation': None,
                'groundtruth_equation': None,
                'R2': None,
                'predicted_complexity': None,
                'groundtruth_complexity': None,
                'initial_conditions': None,
                'ned': None,
                }
            print(f"Error during evaluation: {e}")

        return results








