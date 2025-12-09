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

from utils import utils
from .base import BaseVerifier
from odeformer.model.mixins import FiniteDifferenceMixin
from torchdiffeq import odeint

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
        
        # import ipdb; ipdb.set_trace()
        
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


        return self.candidate_expressions


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


    # def fit_score(self, exprs):
    #     """
    #     Score the expressions based on the number of variables and equations,
    #     with a 60s hard timeout on the fitting step.
    #     """
    #     try:
    #         fit_exprs = self.fit(exprs, self.init_conds, self.time_series)
    #     except MyTimeoutError:
    #         self.logger.info(f"fit_score: fitting timed out after {self.timeout_seconds}s")
    #         return None
    #     except Exception as e:
    #         self.logger.error(f"fit_score: error in fitting: {e}")
    #         return None

    #     try:
    #         pred_data = self.generate([fit_exprs], [self.init_conds])[0]
    #         r2 = utils.r2_score(self.time_series, pred_data)
    #         complexity = sum(e.count_ops() for e in exprs)
    #     except Exception as e:
    #         self.logger.error(f"fit_score: error in generating data: {e}")
    #         return None

    #     return {
    #         "fit_exprs": str(fit_exprs),
    #         "R2": r2,
    #         "Complexity": np.sum(complexity).item()
    #     }

    # def fit_score(self, exprs):
    #     """
    #     Score the expressions based on the number of variables and equations,
    #     with a 60s hard timeout on the entire fit_score execution.
    #     """
    #     import signal
    #     import numpy as np
    #     from odeformer.utils import MyTimeoutError

    #     # Setup a SIGALRM handler for hard timeout
    #     class TimeoutException(Exception): pass
    #     def _handler(signum, frame):
    #         raise TimeoutException()

    #     timeout =  self.timeout_seconds
    #     old_handler = signal.signal(signal.SIGALRM, _handler)
    #     signal.alarm(timeout)
    #     try:
    #         # perform fitting
    #         fit_exprs = self.fit(exprs, self.init_conds, self.time_series)
    #         # generate predictions and compute metrics
    #           # small delay to ensure signal is set up
    #         pred_data = self.generate([fit_exprs], [self.init_conds])[0]
    #         r2 = utils.r2_score(self.time_series, pred_data)
    #         complexity = sum(e.count_ops() for e in exprs)
    #         return {
    #             "fit_exprs": fit_exprs,
    #             "R2": r2,
    #             "Complexity": np.sum(complexity).item()
    #         }
    #     except TimeoutException:
    #         print(f"[DEBUG] fit_score: caught TimeoutException after {timeout}s")
    #         self.logger.info(f"fit_score: timed out after {timeout}s")
    #         return None
    #     except MyTimeoutError:
    #         print(f"[DEBUG] fit_score: caught MyTimeoutError after decorator timeout {timeout}s")
    #         self.logger.info(f"fit_score: timed out after decorator after {timeout}s")
    #         return None
    #     except Exception as e:
    #         print(f"[DEBUG] fit_score: caught Exception: {e}")
    #         self.logger.error(f"fit_score: error: {e}")
    #         return None
    #     finally:
    #         # clear alarm and restore handler
    #         print("[DEBUG] fit_score: clearing alarm and restoring handler")
    #         signal.alarm(0)
    #         signal.signal(signal.SIGALRM, old_handler)

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
            # only treat the symbol 'c' as a learnable coefficient; leave numeric literals intact
            if expr.is_Symbol and str(expr) == 'c':
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



    # def fit(self, raw_exprs, init_conds, gt_data, refine_method="BFGS"):
    #     import numpy as np, sympy
    #     from scipy.integrate import solve_ivp
    #     from scipy.optimize import minimize
    #     # Common setup: replace 'c' symbols and prepare data
    #     exprs_c, c_syms = self.replace_coefficients(raw_exprs)
    #     var_syms = self.dataset.var_symbols
    #     t_eval_arr = np.array(self.solve_config["t_eval"])
    #     # Adjust gt_data: remove time column if present
    #     X_data = gt_data if gt_data.shape[1] == len(var_syms) else gt_data[:, 1:]
    #     deriv = self.approximate_derivative(X_data, t_eval_arr)
    #     fitted = []
    #     idx_c = 0
    #     x0_linear_coeffs=[]
    #     for i, expr in enumerate(exprs_c):
    #         # pick out the coefficients for this equation
    #         local_cs = [c for c in c_syms[idx_c:] if expr.has(c)]
    #         n_local = len(local_cs)
    #         # build design matrix
    #         basis = [expr.coeff(c) for c in local_cs]
    #         fns = [sympy.lambdify(var_syms, b, "numpy") for b in basis]
    #         Phi = np.vstack([fn(*[X_data[:, j] for j in range(X_data.shape[1])]) 
    #                          for fn in fns]).T
    #         y = deriv[:, i]
    #         # filter out any NaNs
    #         mask = (~np.isnan(Phi).any(axis=1)) & (~np.isnan(y))
    #         Phi_clean = Phi[mask]
    #         y_clean = y[mask]
    #         # least squares
    #         cs, *_ = np.linalg.lstsq(Phi_clean, y_clean, rcond=None)
    #         subs = {local_cs[j]: float(cs[j]) for j in range(n_local)}
    #         fitted.append(expr.subs(subs))
    #         idx_c += n_local
    #         x0_linear_coeffs.extend([float(cs[j]) for j in range(n_local)])  # Collect initial guess for BFGS
        
    #     # import ipdb; ipdb.set_trace()   
    #     if not self.refine_coefficients:
    #         return fitted 
        
    #     # Stage 2: Refine coefficients using BFGS optimization
    #     lambdified_exprs_for_bfgs = [sympy.lambdify(var_syms + c_syms, expr, "numpy") for expr in exprs_c]

    #     def objective_function(current_c_values):
    #         def ode_system_rhs(tt, yy):
    #             return np.array([fn(*yy, *current_c_values) for fn in lambdified_exprs_for_bfgs], dtype=float)
            
    #         with np.errstate(invalid='ignore', divide='ignore'): # Suppress warnings during solve_ivp
    #             solution = solve_ivp(ode_system_rhs, 
    #                                  (t_eval_arr[0], t_eval_arr[-1]), 
    #                                  init_conds, 
    #                                  t_eval=t_eval_arr, 
    #                                  **{k:v for k,v in self.solve_config.items() if k not in ("t_span","t_eval")})
            
    #         if not solution.success or solution.y.shape[1] != len(t_eval_arr):
    #             return 1e6 # Penalty for failed or incomplete integration
            
    #         predicted_Y = solution.y.T
    #         return -utils.r2_score(X_data, predicted_Y) # Minimize negative R² (i.e., maximize R²)

    #     # Prepare initial guess for BFGS
    #     current_initial_guess = x0_linear_coeffs
    #     optimization_result = None

    #     if len(c_syms) > 0:
    #         self.logger.info(f"Attempting BFGS with initial guess from linear fit: {current_initial_guess}")
    #         if len(current_initial_guess) != len(c_syms):
    #             self.logger.warning(f"BFGS attempt skipped: initial guess length mismatch. Guess: {len(current_initial_guess)}, Expected: {len(c_syms)}. Using zeros as fallback for this attempt.")
    #             current_initial_guess = np.zeros(len(c_syms))

    #         optimization_result = minimize(objective_function,
    #                                        current_initial_guess,
    #                                        method="BFGS",
    #                                        options={"gtol": 1e-4, "maxiter": 50}) # Reduced maxiter for speed
    #         if optimization_result.success:
    #             self.logger.info(f"BFGS converged successfully.")
    #         else:
    #             self.logger.warning(f"BFGS optimization did not converge: {optimization_result.message}. Using coefficients from linear fit.")
    #     elif len(c_syms) == 0: # No symbolic coefficients to fit
    #         optimization_result = type('MockRes', (), {'success': True, 'x': []})()


    #     final_coefficients = []
    #     if optimization_result is not None and optimization_result.success:
    #         final_coefficients = optimization_result.x
    #     elif len(c_syms) > 0: # BFGS failed or was not run (e.g. initial guess mismatch, though handled)
    #         self.logger.info("Using coefficients from linear fit as BFGS did not succeed or was not applicable.")
    #         final_coefficients = x0_linear_coeffs
    #     # else: final_coefficients remains empty if c_syms is empty, which is correct.


    #     # Substitute final fitted coefficients back into the expressions
    #     fitted_expressions = []
    #     if len(final_coefficients) == len(c_syms):
    #         subs_dict = {c_syms[j]: float(final_coefficients[j]) for j in range(len(c_syms))}
    #         for expr in exprs_c:
    #             fitted_expressions.append(expr.subs(subs_dict))
    #     else: # Fallback if length mismatch (should be rare with current logic)
    #         self.logger.error(f"Final coefficient length ({len(final_coefficients)}) mismatch with c_syms length ({len(c_syms)}). Returning expressions with initial linear fit coeffs or original if that failed.")
    #         # As a safer fallback, try to use x0_linear_coeffs if they match, else exprs_c (partially fitted or as is)
    #         if len(x0_linear_coeffs) == len(c_syms):
    #             subs_dict_fallback = {c_syms[j]: float(x0_linear_coeffs[j]) for j in range(len(c_syms))}
    #             for expr in exprs_c:
    #                 fitted_expressions.append(expr.subs(subs_dict_fallback))
    #         else: # Absolute fallback
    #             fitted_expressions = exprs_c


    #     return fitted_expressions


    # def fit(self, raw_exprs, init_conds, gt_data):
    #     import numpy as np
    #     import sympy

    #     # import ipdb; ipdb.set_trace()

    #     # Common setup: replace 'c' symbols and prepare data
    #     exprs_c, c_syms = self.replace_coefficients(raw_exprs)
    #     var_syms = self.dataset.var_symbols
    #     t_eval_arr = np.array(self.solve_config["t_eval"])
    #     # Adjust gt_data: remove time column if present
    #     X_data = gt_data if gt_data.shape[1] == len(var_syms) else gt_data[:, 1:]
    #     deriv = self.approximate_derivative(X_data, t_eval_arr)
    #     fitted = []
    #     idx_c = 0
    #     x0_linear_coeffs=[]
    #     for i, expr in enumerate(exprs_c):
    #         # pick out the coefficients for this equation
    #         local_cs = [c for c in c_syms[idx_c:] if expr.has(c)]
    #         n_local = len(local_cs)
    #         # build design matrix
    #         basis = [expr.coeff(c) for c in local_cs]
    #         fns = [sympy.lambdify(var_syms, b, "numpy") for b in basis]
    #         Phi = np.vstack([fn(*[X_data[:, j] for j in range(X_data.shape[1])]) 
    #                          for fn in fns]).T
    #         y = deriv[:, i]
    #         # filter out any NaNs
    #         mask = (~np.isnan(Phi).any(axis=1)) & (~np.isnan(y))
    #         Phi_clean = Phi[mask]
    #         y_clean = y[mask]
    #         # least squares
    #         cs, *_ = np.linalg.lstsq(Phi_clean, y_clean, rcond=None)
    #         # cs = [round(float(t), 1) for t in cs]  # Ensure cs is a list of floats
    #         subs = {local_cs[j]: float(cs[j]) for j in range(n_local)}
    #         fitted.append(expr.subs(subs))
    #         idx_c += n_local
    #         x0_linear_coeffs.extend([float(cs[j]) for j in range(n_local)])  # Collect initial guess for BFGS
    #     return fitted
    

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
        # def simplify_node(node):
        #     # If the node is a power operation and the base is a constant in the list
        #     if node.is_Pow and node.base.is_Symbol and node.base not in var_list and node.exp.is_Number:
        #         return node.base  # Replace c**n with c
            
        #     # If the node is a unary function applied to a constant in the list
        #     if node.is_Function and len(node.args) == 1 and node.args[0].is_Symbol and node.args[0] not in var_list:
        #         return  node.args[0]  # Replace func(c) with c
            
        #     if len(node.args)==2 and node.args[0].is_Symbol and node.args[1].is_Symbol and node.args[0] not in var_list and node.args[1] not in var_list:
        #         return  node.args[0]  # Replace func(c, c) with c
        #     # Recursively simplify arguments
        #     if node.is_Atom:
        #         return node
        #     new_args = []
        #     flag=1
        #     for arg in node.args:
        #         simplified_arg = simplify_node(arg)
        #         if not simplified_arg.is_Symbol:
        #             new_args.append(simplified_arg)
        #         elif simplified_arg in var_list:
        #             new_args.append(simplified_arg)
        #         elif simplified_arg not in var_list and flag:
        #             new_args.append(simplified_arg)
        #             flag=0

        #     if len(new_args) == 0:
        #         return sympy.Symbols("c")

        #     return node.func(*new_args) #node.func(*[simplify_node(arg) for arg in node.args])


        def simplify_node(node):
            # Recursively simplify arguments
            if node.is_Atom:
                if node in var_list:
                    return node
                else:
                    return None

            new_args = []
            for arg in node.args:
                if arg.is_Symbol and arg not in var_list:
                    continue
                simplified_arg = simplify_node(arg)
                if simplified_arg is not None:
                    new_args.append(simplified_arg)
                else:
                    continue

            if len(new_args) == 0:
                return None

            return node.func(*new_args) #node.func(*[simplify_node(arg) for arg in node.args])
        if type(expr) is sympy.core.add.Add:
            len_args = len(expr.args)
        else:
            len_args = 1
        # import ipdb; ipdb.set_trace()
        new_expr = simplify_node(expr)

        while new_expr != expr and new_expr is not None:
            expr = new_expr
            new_expr = simplify_node(expr)
        
        local_cs = []
        if new_expr is None:
            new_expr = sympy.symbols("c0")
            return new_expr, [sympy.symbols("c0")]
        idx=0
        if type(new_expr) is sympy.core.add.Add:
            new_args = []
            for i, term in enumerate(new_expr.args):
                c = sympy.symbols(f"c{idx}")
                idx+=1
                local_cs.append(c)
                term = sympy.Mul(c, term)
                new_args.append(term)
            new_expr = sympy.Add(*new_args)
        else:
            new_expr = sympy.Mul(sympy.symbols("c0"), new_expr)
            local_cs.append(sympy.symbols("c0"))
            idx+=1

        current_len = len(new_expr.args) if type(new_expr) is sympy.core.add.Add else 1

        if current_len < len_args:
            local_cs.append(sympy.symbols(f"c{idx}"))
            new_expr = new_expr + sympy.symbols(f"c{idx}")
            

        # import ipdb; ipdb.set_trace()
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



    # def fit(self, raw_exprs, init_conds, gt_data):
    #     import ipdb; ipdb.set_trace()
    #     # fast linear regression on finite-difference derivatives
    #     exprs_c, c_syms = self.replace_coefficients(raw_exprs)
    #     var_syms = self.dataset.var_symbols
    #     # prepare data and time grid
    #     t = np.array(self.solve_config["t_eval"])
    #     X = gt_data if gt_data.shape[1] == len(var_syms) else gt_data[:,1:]
    #     # approximate derivatives once
    #     deriv = self.approximate_derivative(X, t)
    #     fitted = []
    #     # fit each expression independently via least squares
    #     idx_c = 0
    #     for i, expr in enumerate(exprs_c):
    #         # gather local coefficients for this eqn
    #         local_cs = [c for c in c_syms[idx_c:] if expr.has(c)]
    #         n_local = len(local_cs)
    #         # build design matrix Phi of shape (T, n_local)
    #         basis = [expr.coeff(c) for c in local_cs]
    #         # lambdify basis functions (numpy) and solve least squares on GPU via PyTorch
    #         fns = [sympy.lambdify(var_syms, b, 'numpy') for b in basis]
    #         # evaluate design matrix on CPU then move to GPU
    #         Phi_t = torch.cat([
    #             torch.from_numpy(np.atleast_1d(fn(*[X[:, j] for j in range(X.shape[1])]))).to(self.device).unsqueeze(1)
    #             for fn in fns
    #         ], dim=1)
    #         y_t = torch.from_numpy(deriv[:, i]).to(self.device)
    #         # compute least-squares solution via pseudo-inverse on GPU
    #         cs_t = torch.linalg.pinv(Phi_t) @ y_t
    #         cs = cs_t.cpu().numpy()
    #         # substitute back and append
    #         subs = {local_cs[j]: float(cs[j]) for j in range(n_local)}
    #         fitted.append(expr.subs(subs))
    #         idx_c += n_local
    #     import ipdb; ipdb.set_trace()
    #     return fitted



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