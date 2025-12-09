import pandas as pd
import numpy as np
import sympy
import utils
import re
from scipy.integrate import solve_ivp
from tqdm import tqdm
import torch
import json
from copy import copy
import signal
from tqdm import tqdm
from multiprocessing import Process, Queue
import warnings
warnings.filterwarnings("ignore")

from utils import utils
from .base import BaseVerifier
from odeformer.model.mixins import FiniteDifferenceMixin
from torchdiffeq import odeint

from gplearn.genetic import SymbolicRegressor
from gplearn.functions import make_function
from gplearn._program import _Program
from gplearn.functions import _function_map


from models.judge_model import JudgeModel

from odeformer.utils import timeout, MyTimeoutError


class DE_Verifier(BaseVerifier, FiniteDifferenceMixin):
    def __init__(self, cfg, dataset, logger, device=None, dtype=None):
        super().__init__(cfg)
        # initialize finite-difference derivative approximation
        FiniteDifferenceMixin.__init__(self)
        # import ipdb; ipdb.set_trace()
        self.dataset = dataset
        self.logger = logger
        self.device = device
        self.dtype=dtype
        self.judge_model = None
        if cfg.experiment.symbolic_expression.judge.enabled:
            self.judge_name = cfg.experiment.symbolic_expression.judge.get("name", "gpt-4.1-nano")
            # import ipdb; ipdb.set_trace()
            self.judge_model =JudgeModel(self.judge_name, device, dtype)
            judge_prompt_path = cfg.experiment.symbolic_expression.judge.get("prompt_path", None)
            if judge_prompt_path is not None:
                with open(judge_prompt_path, 'r') as f:
                    self.judge_prompt = json.load(f)
                # import ipdb; ipdb.set_trace()
                self.judge_prompt = '\n '.join(self.judge_prompt.get("judge_prompt"))

        # import ipdb; ipdb.set_trace()

        self.use_context = cfg.experiment.symbolic_expression.get("use_context", False)
        self.reasoning = cfg.experiment.symbolic_expression.get("reasoning", False)
        self.additional_prompt_type = cfg.experiment.symbolic_expression.additional_prompt


        self.domain = self.dataset.domain
        self.num_vars = dataset.num_vars
        self.num_exprs = dataset.num_eqs
        self.var_list = dataset.var_symbols
        
        self.gt_exprs = dataset.eqs
        self.init_conds = dataset.initial_conditions
        self.consts_range = dataset.consts_range
    
        self.solve_config = {
            "t_span": (0, 10),
            "method": "LSODA",       # switch to BDF to avoid LSODA prints
            "rtol": 1e-6,
            "atol": 1e-6,
            "first_step": 1e-6,
            "t_eval": np.linspace(0, 10, 150),
            "dense_output": False,
        }
        # no min_step needed for BDF
        if cfg.experiment.symbolic_expression.context.time_series:
            self.time_series = self.generate([self.gt_exprs], [self.init_conds])[0]
        else:
            self.time_series = None
        if cfg.experiment.symbolic_expression.context.var_description:
            self.var_description = dataset.var_description
        else:
            self.var_descripton = None
        if cfg.experiment.symbolic_expression.context.consts_description:
            self.const_description = dataset.const_desciption
        else:
            self.const_description = None

        self.image = dataset.image

        if cfg.get("noise_level", 0) > 0:
            self.time_series = self.inject_noise(self.time_series, cfg.noise_level)

        
        if self.additional_prompt_type == 'time_series' and self.time_series is not None:
            self.additional_prompt = self.time_series
        elif self.additional_prompt_type == 'image' and self.image is not None:
            self.additional_prompt = self.image
        else:
            self.additional_prompt=None
        
        self.tolerance = cfg.experiment.symbolic_expression.get("tolerance", 0.99)

        # import ipdb; ipdb.set_trace()

        self.candidate_expressions = pd.DataFrame(columns=["expr", "R2", "Complexity"])

        self.OOD = cfg.experiment.symbolic_expression.get("OOD", False)

        self.limit_candidates_show=5

        self.timeout_seconds = 15  # adjust as needed
        # self.fit = timeout(seconds=self.timeout_seconds)(self.fit)

        self.max_data_len = 150

        self.patience_limit = 10
        self.patience = self.patience_limit

        # self.scored_expressions = []
    def inject_noise(self, data, noise_level=0.01):
        """
        Add Gaussian noise to time series data.
        
        Args:
            data: Time series data array
            noise_level: Standard deviation of noise as fraction of signal std
        
        Returns:
            Noisy data
        """
        # Scale noise relative to signal standard deviation for each variable
        signal_std = np.std(data, axis=0, keepdims=True)
        noise = noise_level * signal_std * np.random.randn(*data.shape)
        return data + noise


    def prepare_prompt(self, base_prompt):
        """
        Prepare the prompt for the model.
        """
        # import ipdb; ipdb.set_trace()
        prompt = '\n '.join(base_prompt)
        try:
            var_points = ""
            var_list = [str(x).replace("_","") for x in self.var_list]
            for i, var in enumerate(var_list):
                time_series = str(np.round(self.time_series[:self.max_data_len,i],3).tolist())
                if self.use_context and self.additional_prompt_type == 'time_series' and self.time_series is not None:
                    var_points += f'{var}: <ts><ts/>; '
                    # import ipdb; ipdb.set_trace()
                elif self.use_context and self.additional_prompt_type == 'image' and self.image is not None:
                    var_points += f"{var}: variable {i+1} in the image; "
                elif self.use_context:
                    var_points += f"{var}: {time_series};"
            self.var_points = str(var_points)
            prompt = prompt.replace("{var_points}", str(var_points))
            prompt = prompt.replace("{num_vars}", str(len(var_list)))
            prompt = prompt.replace("{num_eqs}", str(self.num_exprs))
            prompt = prompt.replace("{var_list}", str(var_list))

        except:
            print("No time series data input.")

        if self.use_context:
            context = ""
            if self.domain is not None:
                context += f"The domain of this problem is called '{self.domain}'.\n"
            context += f"Here are the variables and their descriptions: {self.var_description}.\n"
            self.context = context

            prompt = prompt.replace("{context}", context)


        if not self.candidate_expressions.empty:
            prompt = prompt.replace("{functions}", str(self.candidate_expressions.iloc[:self.limit_candidates_show]))

        if self.num_exprs == 1:
            example_str = r'{"eq": "c*x_0 + c",  "dim": 1'
        elif self.num_exprs == 2:
            example_str = r'{"eq": "c*x_0 + c*x_1 + c | c*x_0 - c*x_1",  "dim": 2'
        elif self.num_exprs == 3:
            example_str = r'{"eq": "c*x_0 - c*x_1 + c*x_2 + c | c*x_1*x_2 - c*x_0 | c*x_0 + c*x_1 - c*x_2",  "dim": 3'
        elif self.num_exprs == 4:
            example_str = r'{"eq": "c*x_0*x_1 + c*x_2 - c*x_3 + c | c*x_1/x_0 + c*x_3 | c*x_2 - c*x_0*x_3 | c*x_3 + c*x_0 - c*x_1",  "dim": 4'
        
        if self.reasoning:
            example_str += r', "reasoning": "...(explain why the generated equations fits)"}'
        else:
            example_str += r'}'
        # import ipdb; ipdb.set_trace()
        prompt = prompt.replace("{example_str}", example_str)

        self.current_prompt = prompt

        return prompt

        
    
    def generate(self, exprs_list, init_conds_list):
        """
        this code is adopeted from odebench
        """
        # import ipdb; ipdb.set_trace()
        solutions = []
        for i, exprs in tqdm(enumerate(exprs_list)):
            var_symbols = self.dataset.var_symbols
            # import ipdb; ipdb.set_trace()
            callable_fn = lambda t, x: np.array([f(*x) for f in [sympy.lambdify(var_symbols, eq, 'numpy') for eq in exprs]])
            init_conds = init_conds_list[i]
            # suppress invalid-value-in-power warnings during integration
            with np.errstate(invalid='ignore', divide='ignore'):
                # import ipdb; ipdb.set_trace()
                sol = solve_ivp(callable_fn, **self.solve_config, y0=init_conds)
            # import ipdb; ipdb.set_trace()
            sol_dict = {
                "success": sol.success,
                "message": sol.message,
                "t": sol.t.tolist(),
                "y": sol.y.tolist(),
                "nfev": int(sol.nfev),
                "njev": int(sol.njev),
                "nlu": int(sol.nlu),
                "status": int(sol.status),
            }
            if sol.status != 0:
                print(f"Error in expression set {i}: {sol.message}")
            solutions.append(np.array(sol_dict['y']).T)
        
        return solutions
        
    def parse_model_output(self, str_funcs):
        """
        Parse model output in JSON-like list of dicts format, extracting and cleaning each equation.
        Maintains "|" separators if present in the equations.
        Example input:
        'Here are five samples of synthetic ODEs in the required format:\n\n```\n{\n   "eq": "c*x0|c*x1",\n   "dim": 2\n}, ...'
        Returns a list of lists of sympy expressions (one per sample).
        """
        # import ipdb; ipdb.set_trace()
        # expr = 'c*x_0  + c * x_0 * x_1 | c*x_1 + c * x_0*x_1'
        # expr = self.dataset.str_to_sympy(expr, var_list=self.var_list)
        # fitted = self.fit(expr, self.init_conds, self.time_series)
        # score = utils.r2_score(self.time_series, self.generate([fitted], [self.init_conds])[0])
        # import ipdb; ipdb.set_trace()

        import re

        # Extract all dict-like blocks
        dict_blocks = re.findall(r'\{[^}]*\}', str_funcs)
        all_expressions = []
        reasonings = []

        for block in dict_blocks:
            # Extract the equation string after "eq":
            # eq_match = re.search(r'"eq"\s*:\s*"([^"]+)"', block)
            eq_match = re.search(r"""['"]eq['"]\s*:\s*['"]([^'"]+)['"]""", block)
            if not eq_match:
                self.logger.info(f"Skipping block as no 'eq' found: {block}")
                continue
            eq_str = eq_match.group(1)
            # Split by "|" but keep the separator in the string for later joining
            eqs = [e.strip() for e in eq_str.split("|")]
            cleaned_eqs = []
            for eq in eqs:
                cleaned = self.clean_function(eq)
                if cleaned != "":
                    cleaned_eqs.append(cleaned)
            # Re-join with "|" if there were multiple equations
            joined_cleaned = " | ".join(cleaned_eqs)
            if joined_cleaned:
                try:
                    # If you want to keep as string with "|", append as is
                    # import ipdb; ipdb.set_trace()
                    function = self.dataset.str_to_sympy(joined_cleaned, var_list=self.var_list)
                    # import ipdb; ipdb.set_trace()
                    if len(function) != self.num_exprs:
                        self.logger.info(f"Skipping block as number of equations does not match: {block}")
                        continue
                    all_expressions.append(function)
                except Exception as e:
                    self.logger.warning(f"Could not parse line {joined_cleaned}.")
                    self.logger.warning(str(e))
                    continue
            
            reasoning_match = re.search(r'"reasoning"\s*:\s*"([^"]*)"', block)
            # import ipdb; ipdb.set_trace()
            reasonings.append(reasoning_match.group(1) if reasoning_match else None)


        return all_expressions, reasonings
    

    def get_best_expr(self,):
        return self.candidate_expressions.iloc[0]['expr'], self.candidate_expressions.iloc[0]['fit_exprs'], self.candidate_expressions.iloc[0]['R2']

    def check_tolerance(self,):
        """
        Check if the candidate expressions are within the tolerance.
        """
        if self.candidate_expressions.empty:
            return False
        # import ipdb; ipdb.set_trace()

        if self.candidate_expressions.iloc[0]["R2"] >= self.tolerance:
            return True
        else:
            return False
        


    def check_duplicate(self, expr):
        """
        Check if the expression is already in the candidate expressions.
        """
        # if self.candidate_expressions.empty:
        #     return False
        # # import ipdb; ipdb.set_trace()
        # for candidate in self.candidate_expressions:
        #     if expr == candidate.expr.item():
        #         return True
        return False
    
    def update_pareto_front(self,):
        if self.candidate_expressions.empty:
            return
        # import ipdb; ipdb.set_trace()
        if 'context_alignment' in self.candidate_expressions.columns:
            self.candidate_expressions.sort_values(by=["R2", "context_alignment", "Complexity"], ascending=[False, False, True], inplace=True)
        else:
            self.candidate_expressions.sort_values(by=["R2", "Complexity"], ascending=[False, True], inplace=True)
        self.candidate_expressions.drop_duplicates(subset=["expr"], keep='first', inplace=True)

        pareto_indices_loc = [] 
        for i in range(len(self.candidate_expressions)):
            candidate_row = self.candidate_expressions.iloc[i]
            is_dominated = False
            
            for p_loc in pareto_indices_loc:
                pareto_member_row = self.candidate_expressions.iloc[p_loc]
                
                if (pareto_member_row["R2"] >= candidate_row["R2"] and \
                    pareto_member_row["Complexity"] <= candidate_row["Complexity"] and \
                    pareto_member_row['context_alignment'] >= candidate_row['context_alignment'] if 'context_alignment' in self.candidate_expressions.columns else True
                    ):
                    is_dominated = True
                    break
            
            if not is_dominated:
                pareto_indices_loc.append(i)
        
        if frozenset(self.candidate_expressions['expr']) == frozenset(self.candidate_expressions.iloc[pareto_indices_loc]['expr']):
            self.patience -= 1
        else:
            self.patience = self.patience_limit
        
        self.candidate_expressions = self.candidate_expressions.iloc[pareto_indices_loc].reset_index(drop=True)


    def update_candidates(self, expr_list, hybrid=False):
        """
        Update the candidate expressions with the new expression.
        """
        # import ipdb; ipdb.set_trace()
        for expr, reasoning in tqdm(tuple(expr_list)):
            if not self.check_duplicate(expr):
                # import ipdb; ipcdb.set_trace()
                fit_score = self.fit_score(expr)
                if fit_score is None:
                    self.logger.info(f"Expression {expr} is not valid.")
                    continue
                # import ipdb; ipdb.set_trace()
                tmp  ={"expr": str(expr), "reasoning": reasoning}
                tmp.update(**fit_score)
                if self.judge_model is not None:
                    tmp = self.score_quality(tmp)
                tmp = pd.DataFrame({k: [v] for k, v in tmp.items()})
                self.candidate_expressions = pd.concat([self.candidate_expressions, tmp], ignore_index=True)
            else:
                self.logger.info(f"Expression {expr} already exists in the candidate expressions.")

        # self.update_pareto_front() 
        if 'context_alignment' in self.candidate_expressions.columns:
            self.candidate_expressions.sort_values(by=["R2", "context_alignment", "Complexity"], ascending=[False, False, True], inplace=True)
        else:
            self.candidate_expressions.sort_values(by=["R2", "Complexity"], ascending=[False, True], inplace=True)
        
        # import ipdb; ipdb.set_trace()        
        self.candidate_expressions.drop_duplicates(subset=["expr"], keep='first', inplace=True)
        self.candidate_expressions = self.candidate_expressions.head(100)
        self.logger.info(f"Num candidate expressions: {len(self.candidate_expressions)}")

        if hybrid:
            self.expand_candidates()

        # import ipdb; ipdb.set_trace()

        return self.candidate_expressions


    def gplearn_to_sympy(self, expr_str):
        # collect all X-variable names (X0, X1, …)
        from sympy import symbols, sympify, sin, cos, log
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


    def _sympy_to_gplearn(self, expr, est_gp):
        """
        Converts a SymPy expression to a gplearn program list by traversing the expression tree.
        This version correctly handles n-ary SymPy functions by converting them to nested binary gplearn functions.
        """
        # Map function names to the actual _Function objects from the fitted estimator
        function_map = {f: _function_map[f] for f in est_gp.function_set}
        feature_map = {name: i for i, name in enumerate(est_gp.feature_names)}

        # This mapping helps translate SymPy function classes to gplearn function names
        sympy_to_gplearn_map = {
            sympy.Add: 'add',
            sympy.Mul: 'mul',
            sympy.Pow: 'pow',
            sympy.sin: 'sin',
            sympy.cos: 'cos',
            sympy.tan: 'tan',
            sympy.Abs: 'abs',
            sympy.log: 'log',
            sympy.sqrt: 'sqrt',
            sympy.Min: 'min',
            sympy.Max: 'max',
            # Add core types for robustness
            sympy.core.add.Add: 'add',
            sympy.core.mul.Mul: 'mul',
            sympy.core.power.Pow: 'pow',
        }

        program = []
        constants = []

        def traverse(sub_expr):
            """
            Recursively traverses the SymPy expression tree (depth-first) and builds the gplearn program list.
            """
            # Case 1: The node is a function or an operator (e.g., Add, Mul, sin)
            if not sub_expr.is_Atom:
                gplearn_name = sympy_to_gplearn_map.get(type(sub_expr))
                if not gplearn_name or gplearn_name not in function_map:
                    raise ValueError(f"Unsupported SymPy function: {type(sub_expr)}")

                gplearn_func = function_map[gplearn_name]
                args = sub_expr.args

                # Handle n-ary Add/Mul by converting to a nested binary tree
                if gplearn_name in ('add', 'mul') and len(args) > 2:
                    # Convert Add(a, b, c) to add(a, add(b, c))
                    program.append(gplearn_func)
                    traverse(args[0])
                    # Create a new SymPy expression for the rest of the arguments
                    remaining_expr = sub_expr.func(*args[1:])
                    traverse(remaining_expr)
                else:
                    # Handle standard binary and unary functions
                    if len(args) != gplearn_func.arity:
                        raise ValueError(f"Arity mismatch for {gplearn_name}: expected {gplearn_func.arity}, got {len(args)}")
                    program.append(gplearn_func)
                    for arg in args:
                        traverse(arg)
                return

            # Case 2: The node is a variable (e.g., x_0, x_1)
            if sub_expr.is_Symbol:
                if str(sub_expr) in feature_map:
                    program.append(feature_map[str(sub_expr)])
                else:  # Treat as a learnable constant
                    program.append(-1)
                    constants.append(1.0)
                return

            # Case 3: The node is a numeric literal (e.g., 3.14, -2)
            if sub_expr.is_number:
                program.append(-1)
                constants.append(float(sub_expr))
                return

            raise ValueError(f"Unsupported SymPy expression type: {type(sub_expr)}")

        # Start the traversal from the root of the expression tree
        # import ipdb; ipdb.set_trace()
        traverse(expr)
        return program, constants
    
    def validate_program(self, program):
        from gplearn.functions import _Function
        """Rough check that the embedded program in the object is valid."""
        terminals = [0]
        for node in program:
            if isinstance(node, _Function):
                terminals.append(node.arity)
            else:
                terminals[-1] -= 1
                while terminals[-1] == 0:
                    terminals.pop()
                    terminals[-1] -= 1
        return terminals == [-1]

    def expand_candidates(self):
        # use the top N candidates as initial popluation and use Genetic programming to expand the candidates
        
        self.logger.info("Expanding candidates using Genetic Programming.")
        top_candidates = self.candidate_expressions.head(5)
        if top_candidates.empty:
            self.logger.info("No candidates to expand.")
            return

        # Prepare data for GP
        t_eval_arr = np.array(self.solve_config["t_eval"])
        X_data = self.time_series if self.time_series.shape[1] == len(self.var_list) else self.time_series[:, 1:]
        deriv = self.approximate_derivative(X_data, t_eval_arr)

        new_expressions = []
        
        # Define the symbolic regressor template
        est_gp_template = SymbolicRegressor(population_size=len(top_candidates),
                                   generations=5, stopping_criteria=0.01,
                                   p_crossover=0.7, p_subtree_mutation=0.1,
                                   p_hoist_mutation=0.05, p_point_mutation=0.1,
                                   max_samples=0.9, verbose=0,
                                   parsimony_coefficient=0.01, random_state=0,
                                   feature_names=[str(v) for v in self.var_list],
                                   function_set=('add', 'sub', 'mul', 'div', 'sin', 'cos', 'log'))
        # Fit with dummy data to initialize internal structures like function_set_
        est_gp_template.fit(np.zeros((2, len(self.var_list))), np.zeros(2))


        for i in range(self.num_exprs):
            y_target = deriv[:, i]
            # Use top candidates for this equation as initial programs
            initial_programs = []
            for _, row in top_candidates.iterrows():
                try:
                    # The expression is stored as a string, so we need to parse it back to sympy
                    expr_obj = sympy.sympify(row['fit_exprs'][i])
                    
                    # Convert sympy expression to gplearn program list
                    program_list, constants = self._sympy_to_gplearn(expr_obj, est_gp_template)
                    
                    # Create a _Program object
                    program = _Program(function_set=est_gp_template.function_set,
                                       arities=est_gp_template._arities,
                                       init_depth=est_gp_template.init_depth, 
                                       init_method=None, # Not needed when program is provided
                                       n_features=est_gp_template.n_features_in_,
                                       const_range=est_gp_template.const_range,
                                       metric=est_gp_template._metric,
                                       p_point_replace=est_gp_template.p_point_replace,
                                       parsimony_coefficient=est_gp_template.parsimony_coefficient,
                                       random_state=None, # Let the regressor manage state
                                       feature_names=est_gp_template.feature_names,
                                       program=program_list)
                    program.constants = np.array(constants)
                    # import ipdb; ipdb.set_trace()
                    y_pred = program.execute(X_data)
                    program._raw_fitness = est_gp_template._metric(y_target, y_pred, None, None)
                    program.fitness_ = program._raw_fitness + (program.parsimony_coefficient * program.length_)

                    initial_programs.append(program)
                except Exception as e:
                    self.logger.warning(f"Could not parse or convert expression for GP init: {row['expr']}. Error: {e}")

            # Create a new regressor for the actual fitting
            est_gp = SymbolicRegressor(population_size=len(top_candidates) if initial_programs else 50,
                                       generations=50, stopping_criteria=0.01,
                                       p_crossover=0.7, p_subtree_mutation=0.1,
                                       p_hoist_mutation=0.05, p_point_mutation=0.1,
                                       max_samples=0.9, verbose=0,
                                       parsimony_coefficient=0.01, random_state=0,
                                       feature_names=[str(v) for v in self.var_list],
                                       warm_start=False,
                                       function_set=('add', 'sub', 'mul', 'div', 'sin', 'cos', 'log'))

            if not initial_programs:
                self.logger.warning(f"No valid initial programs for equation {i}, using random population.")
                est_gp.fit(X_data, y_target)
            else:
                self.logger.info(f"Running GP for equation {i+1}/{self.num_exprs} with {len(initial_programs)} initial programs...")
                est_gp._programs = [initial_programs]
                est_gp.fit(X_data, y_target)

            try:
                # The best program is a string expression
                new_expr_str = str(est_gp._program)
                new_expressions.append(new_expr_str)
            except Exception as e:
                self.logger.error(f"GP failed for equation {i}: {e}")
                # As a fallback, reuse the best existing expression for this dimension
                new_expressions.append(top_candidates.iloc[0]['fit_exprs'][i])
        
        # import ipdb; ipdb.set_trace()
        # add code here to do the final evaluation on R2 score: pred = self.generate([fit_exprs], [init_conds])[0]; r2 = utils.r2_score(self.time_series, pred)
        if new_expressions:
            try:
                # Convert gplearn string expressions to sympy expressions
                new_sympy_exprs = [self.gplearn_to_sympy(expr_str) if isinstance(expr_str, str) else expr_str for expr_str in new_expressions]
                self.timeout_seconds = 60
                # Score the new candidate expression
                fit_score = self.fit_score(new_sympy_exprs)
                self.timeout_seconds = 15


                if fit_score:
                    # Create a new candidate entry
                    new_candidate = {
                        "expr": " | ".join(map(str, new_sympy_exprs)),
                        "reasoning": "Generated by Genetic Programming",
                        **fit_score
                    }
                    if self.judge_model is not None:
                        new_candidate = self.score_quality(new_candidate)

                    # Add the new candidate to the DataFrame
                    new_candidate_df = pd.DataFrame({k: [v] for k, v in new_candidate.items()})
                    self.candidate_expressions = pd.concat([self.candidate_expressions, new_candidate_df], ignore_index=True)
                    self.logger.info(f"Added new candidate from GP with R2: {fit_score.get('R2')}")
                else:
                    self.logger.warning("Scoring of GP-generated expression failed.")

            except Exception as e:
                self.logger.error(f"Error processing GP-generated expressions: {e}")



    def score_quality(self, scored_expr):
        # Prepare and send prompt to judge model
        base_prompt = copy(self.judge_prompt)
        filled_prompt = base_prompt.replace("{data}", self.var_points)
        filled_prompt = filled_prompt.replace("{context}", self.context)
        filled_prompt = filled_prompt.replace("{candidate_exprs}", str(scored_expr))
        # import ipdb; ipdb.set_trace()
        str_scores = self.judge_model.generate(filled_prompt)

        # Initialize fixed result structure with required keys
        if self.reasoning:
            final = {
                "context_alignment": 0,
                "scientific_plausibility": 0,
                "conciseness_and_clarity": None,
                "logical_coherence": None
            }
        else:
            final = {
                "context_alignment": 0,
                "scientific_plausibility": 0,
            }
        # Parse numeric scores into raw dict (try JSON then regex)
        raw = {}
        try:
            jm = re.search(r'(\{.*?\})', str_scores, re.DOTALL)
            if jm:
                for k, v in json.loads(jm.group(1)).items():
                    iv = int(v)
                    if not (1 <= iv <= 5): raise ValueError
                    raw[k] = iv
            else:
                raise ValueError
        except Exception:
            for key, val_str in re.findall(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*:\s*([0-9]+)', str_scores):
                iv = int(val_str)
                if 1 <= iv <= 5 and key not in raw:
                    raw[key] = iv
        # Map only the predefined keys
        for k in final.keys():
            if k in raw:
                final[k] = raw[k]
        # Merge into scored_expr and return
        # import ipdb; ipdb.set_trace()
        scored_expr.update(final)
        return scored_expr


    def fit_score(self, exprs):
        """
        Score the expressions by running fit+generate in a subprocess, terminate after timeout.
        """

        def worker(exprs, init_conds, time_series, out_q):
            try:
                fit_exprs = self.fit(exprs, init_conds, time_series)
                pred = self.generate([fit_exprs], [init_conds])[0]
                r2 = utils.r2_score(self.time_series, pred)
                comp = sum(e.count_ops() for e in exprs)
                out_q.put(("ok", (fit_exprs, r2, float(np.sum(comp)))))
            except Exception as e:
                out_q.put(("err", e))

        q = Queue()
        p = Process(target=worker, args=(exprs, self.init_conds, self.time_series, q))
        p.start()
        p.join(self.timeout_seconds)
        if p.is_alive():
            p.terminate()
            self.logger.info(f"fit_score: timed out after {self.timeout_seconds}s")
            return None

        if q.empty():
            self.logger.error("fit_score: no result returned")
            return None

        status, payload = q.get()
        if status == "ok":
            fit_exprs, r2, comp = payload
            return {"fit_exprs": fit_exprs, "R2": r2, "Complexity": comp}
        else:
            self.logger.error(f"fit_score: error: {payload}")
            return None
    


    def replace_coefficients(self, exp):
        """
        Replace every numeric constant and every occurrence of 'c' in the expression(s)
        with a unique symbol (c0, c1, ...), even if they appear multiple times.
        Returns (new_expressions, c_symbols) where c_symbols is the list of all new sympy symbols in order of appearance.
        """
        import sympy
        if not isinstance(exp, list):
            exp = [exp]
        counter = [0]
        c_syms = []
        def repl(expr):
            if expr.is_number or (expr.is_Symbol and str(expr) == 'c'): # 
                sym = sympy.Symbol(f'c{counter[0]}')
                c_syms.append(sym)
                counter[0] += 1
                return sym
            elif expr.is_Atom:
                return expr
            else:
                return expr.func(*[repl(arg) for arg in expr.args])
        new_exp = [repl(e) for e in exp]
        return new_exp, c_syms

    

    def fit(self, raw_exprs, init_conds, gt_data):
        import numpy as np
        import sympy

        # Common setup: replace 'c' symbols and prepare data
        exprs_c, c_syms = self.replace_coefficients(raw_exprs)
        var_syms = self.dataset.var_symbols
        t_eval_arr = np.array(self.solve_config["t_eval"])
        # Adjust gt_data: remove time column if present
        X_data = gt_data if gt_data.shape[1] == len(var_syms) else gt_data[:, 1:]
        deriv = self.approximate_derivative(X_data, t_eval_arr)
        fitted = []
        idx_c = 0
        x0_linear_coeffs=[]

        # import ipdb; ipdb.set_trace()
        for i, expr in enumerate(exprs_c):
            # pick out the coefficients for this equation
            # import ipdb; ipdb.set_trace()
            expr, local_cs = self.simplify_constants(expr.expand(), var_syms)
            # expr = self.move_constants_to_front(expr, var_syms)
            # tmp = sympy.parse_expr("c0*x_0 + c1*x_1 + c2*sin(x_2) + c3")
            # local_cs = [c for c in c_syms[idx_c:] if expr.has(c)]
            n_local = len(local_cs)

            basis = [expr.coeff(c) for c in local_cs]
            # basis = self.extract_basis(expr, local_cs)
            fns = [sympy.lambdify(var_syms, b, "numpy") for b in basis]
            
            # Phi = np.vstack([fn(*[X_data[:, j] for j in range(X_data.shape[1])]) 
            #                  for fn in fns]).T
            Phi=[]
            for fn in fns:
                # Check if the function returns a scalar or an array
                result = fn(*[X_data[:, j] for j in range(X_data.shape[1])])
                if not np.isscalar(result):
                    Phi.append(result)
                else:
                    Phi.append(np.full(X_data.shape[0], result))
            Phi = np.vstack(Phi).T

            y = deriv[:, i]
            # least squares
            # import ipdb; ipdb.set_trace()
            cs, *_ = np.linalg.lstsq(Phi, y, rcond=None)
            cs = np.round(cs, decimals=4)
            # cs = [round(float(t), 1) for t in cs]  # Ensure cs is a list of floats
            subs = {local_cs[j]: float(cs[j]) for j in range(n_local)}
            fitted.append(expr.subs(subs))
            idx_c += n_local
            x0_linear_coeffs.extend([float(cs[j]) for j in range(n_local)]) 

        return fitted

    def simplify_constants(self, expr, var_list):
        def simplify_node(node):
            # Recursively simplify arguments
            if node.is_Atom:
                if node in var_list:
                    return node
                else:
                    return None

            new_args = []
            for arg in node.args:
                simplified_arg = simplify_node(arg)
                if simplified_arg is not None:
                    new_args.append(simplified_arg)

            # If no arguments remain after simplification, this branch is pruned
            if not new_args:
                return None
            
            # For Pow, if only one argument remains, it's no longer a power operation.
            if isinstance(node, sympy.Pow) and len(new_args) == 1:
                return new_args[0]

            # Reconstruct the node with simplified arguments
            return node.func(*new_args)

        if type(expr) is sympy.core.add.Add:
            len_args = len(expr.args)
        else:
            len_args = 1
        
        simplified_expr = simplify_node(expr)

        local_cs = []
        if simplified_expr is None:
            # If the entire expression simplifies to None, it means it was purely constant.
            # We represent this with a single new constant.
            c = sympy.Symbol('c0')
            local_cs.append(c)
            return c, local_cs

        idx = 0
        if isinstance(simplified_expr, sympy.Add):
            new_args = []
            for term in simplified_expr.args:
                c = sympy.Symbol(f"c{idx}")
                idx += 1
                local_cs.append(c)
                new_args.append(sympy.Mul(c, term))
            new_expr = sympy.Add(*new_args)
        else:
            c = sympy.Symbol('c0')
            local_cs.append(c)
            new_expr = sympy.Mul(c, simplified_expr)
            idx += 1
        
        # Add a constant term if the number of terms was reduced
        current_len = len(new_expr.args) if isinstance(new_expr, sympy.Add) else 1
        if current_len < len_args:
            c = sympy.Symbol(f"c{idx}")
            local_cs.append(c)
            new_expr = new_expr + c

        return new_expr, local_cs


    def extract_basis(self, expr, symbols):
        """
        Extracts terms involving the given symbols from the expression.
        Returns a list of simplified basis functions.
        """
        basis = []
        for symbol in symbols:
            # Extract terms from the expression
            terms = expr.as_ordered_terms() 
            symbol_terms = [term for term in terms if term.has(symbol)]
            # Combine terms into a single basis function
            basis_function = sum(symbol_terms)
            basis.append(basis_function)
        return basis

    def clean_function(self, function: str) -> str:
        """
        Cleans a single function string to be evaluable and processable by sympy, robustly handling LLM output artifacts.
        - Ensures all function calls are in the form func(...)
        - Ensures all variable symbols are from self.var_list
        - Replaces all numeric values with 'c' (at the very end)
        - Ensures all operators are correctly connected
        - Handles 'c' as a constant only at the very last step
        """
        import re
        allowed_functions = [
            'sqrt', 'exp', 'log', 'abs',
            'sin', 'cos', 'tan', 'cot', 'sinh', 'cosh', 'tanh'
        ]
        allowed_vars = set(str(v) for v in self.var_list)
        eq = function.strip()
        # import ipdb; ipdb.set_trace()
        # Remove LaTeX and unicode artifacts, but do NOT remove any brackets
        eq = re.sub(r'\\?frac\{d[ ,]*([a-zA-Z0-9_]+)\}{dt\}?', r'd\1/dt', eq)
        eq = re.sub(r'\\?cdot', '*', eq)
        eq = re.sub(r'\\?mathrm\{([a-zA-Z0-9_]+)\}', r'\1', eq)
        eq = re.sub(r'\\?left|\\?right', '', eq)
        eq = re.sub(r'\\?begin\{.*?\}|\\?end\{.*?\}', '', eq)
        eq = re.sub(r'\\?displaystyle', '', eq)
        eq = eq.replace(' ', '').replace('\u200b', '')
        eq = re.sub(r'x[,\\]+x', 'x0', eq)
        eq = re.sub(r'\\[a-zA-Z]+', '', eq)
        eq = re.sub(r'^=+', '', eq)
        eq = re.sub(r'^:+', '', eq)
        eq = re.sub(r'^[+*/\-]+', '', eq)
        eq = re.sub(r'[+*/\-]+$', '', eq)
        # Fix unclosed brackets: add missing closing brackets for each function if needed
        # Count open and close brackets
        open_count = eq.count('(')
        close_count = eq.count(')')
        if open_count > close_count:
            eq += ')' * (open_count - close_count)
        # Optionally, remove excess closing brackets at the end
        elif close_count > open_count:
            # Only trim trailing unmatched closing brackets
            diff = close_count - open_count
            for _ in range(diff):
                if eq.endswith(')'):
                    eq = eq[:-1]
                else:
                    break
        # Extract right-hand side if in the form x0 = ...
        match = re.match(r'([a-zA-Z_][a-zA-Z0-9_]*)=([^=]+)', eq)
        if match:
            eq = match.group(2)
        elif '=' in eq:
            parts = eq.split('=')
            if len(parts) > 1:
                eq = parts[1]
        # import ipdb; ipdb.set_trace()
        # 1. Standardize variable names (e.g., x0 to x_0)
        def fix_var(match_obj):
            var = match_obj.group(0)
            if var in allowed_vars:
                return var
            m = re.match(r'x([0-9]+)$', var)
            if m:
                candidate = f"x_{m.group(1)}"
                if candidate in allowed_vars:
                    return candidate
            return var
        eq = re.sub(r'[a-zA-Z_][a-zA-Z0-9_]*', fix_var, eq)
        # 2. Replace ^ with ** for sympy
        eq = eq.replace('^', '**')
        # 3. Ensure all allowed functions are followed by brackets (func(...))
        #    - If func is followed by * or directly by a variable/parenthesis, wrap the next token in brackets
        # import ipdb; ipdb.set_trace()
        for func in allowed_functions:
            # func*... or func... (not already func(...))
            eq = re.sub(rf'\b{func}\*\*', rf'{func}**', eq)  # don't break powers
            eq = re.sub(rf'\b{func}\*\s*([a-zA-Z_][a-zA-Z0-9_]*)', rf'{func}(\1)', eq)
            eq = re.sub(rf'\b{func}\*\s*\(([^)]+)\)', rf'{func}(\1)', eq)
            eq = re.sub(rf'\b{func}\s+([a-zA-Z_][a-zA-Z0-9_]*)', rf'{func}(\1)', eq)
            # If func is directly followed by a variable (no *), e.g. sinx_0 -> sin(x_0)
            eq = re.sub(rf'\b{func}(x_\d+)', rf'{func}(\1)', eq)
        # 4. Remove any accidental double multiplication
        # import ipdb; ipdb.set_trace()
        eq = re.sub(r'\*\*+', '**', eq)
        eq = re.sub(r'\*+', '*', eq)
        eq = re.sub(r'\*\)', ')', eq)
        eq = re.sub(r'\(\*', '(', eq)
        eq = re.sub(r'^\*+', '', eq)
        eq = re.sub(r'\*+$', '', eq)
        # 5. Replace all numeric values (integers, floats, scientific notation) with a unique placeholder
        # import ipdb; ipdb.set_trace()
        eq = re.sub(r'(?<![a-zA-Z_])([0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?)', '__CONST__', eq)
        # 6. Only now, insert multiplication for __CONST__ followed by variable or function
        eq = re.sub(r'(__CONST__)([a-zA-Z_][a-zA-Z0-9_\(])', r'\1*\2', eq)
        # 7. Finally, replace all __CONST__ with 'c'
        eq = eq.replace('__CONST__', 'c')
        eq = eq.replace('c**2', 'c')
        return eq