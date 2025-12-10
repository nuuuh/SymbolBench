import pandas as pd
import numpy as np
import sympy
import utils
import re
from scipy.integrate import solve_ivp
from tqdm import tqdm
import torch

from utils import utils
from .base import BaseVerifier
from torchdiffeq import odeint


class DE_Verifier(BaseVerifier):
    def __init__(self, cfg, dataset, logger):
        super().__init__(cfg)
        # import ipdb; ipdb.set_trace()
        self.dataset = dataset
        self.logger = logger

        self.num_vars = dataset.num_vars
        self.num_exprs = dataset.num_eqs
        self.var_list = dataset.var_symbols
        
        self.gt_exprs = dataset.eqs
        self.init_conds = dataset.initial_conditions
        self.consts_range = dataset.consts_range
    
        self.solve_config = {
            "t_span": (0, 10),
            "method": "BDF",       # switch to BDF to avoid LSODA prints
            "rtol": 1e-3,
            "atol": 1e-6,
            "first_step": 1e-6,
            "t_eval": np.linspace(0, 10, 150),
            "dense_output": False,
        }
        # no min_step needed for BDF
        if cfg.experiment.symbolic_expression.context.time_series:
            self.time_series = self.generate([self.gt_exprs], [self.init_conds])[0]
        if cfg.experiment.symbolic_expression.context.var_description:
            self.var_description = dataset.var_description
        if cfg.experiment.symbolic_expression.context.consts_description:
            self.const_description = dataset.const_desciption
        if cfg.experiment.symbolic_expression.context.image:
            self.image = dataset.image
        
        self.tolerance = cfg.experiment.symbolic_expression.get("tolerance", 0.9)

        # import ipdb; ipdb.set_trace()

        self.candidate_expressions = pd.DataFrame(columns=["expr", "R2", "Complexity"])

        self.OOD = cfg.experiment.symbolic_expression.get("OOD", False)

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
                time_series = str(np.round(self.time_series[:,0],3).tolist())
                var_points += f"{var}: {time_series}; "
            prompt = prompt.replace("{var_points}", str(var_points))
            prompt = prompt.replace("{num_vars}", str(len(var_list)))
            prompt = prompt.replace("{num_eqs}", str(self.num_exprs))
            prompt = prompt.replace("{var_list}", str(var_list))

            incomplete_exprs = ""
            for i in range(len(self.gt_exprs)):
                incomplete_exprs += f"dx{i}/dt = x{i} "
                if i != len(self.gt_exprs) - 1:
                    incomplete_exprs += "|"
            prompt = prompt.replace("{incomplete_exprs}", incomplete_exprs)

        except:
            print("No time series data input.")
        if "{var_description}" in prompt:
            prompt = prompt.replace("{var_description}", str(self.var_description))
        else:
            print("No variable description input.")
        if "{const_description}" in prompt:
            prompt = prompt.replace("{const_description}", str(self.const_description))
        else:
            print("No constant description input.")
        # if "{image}" in prompt:
        #     prompt = prompt.replace("{image}", str(self.image))
        # else:
        #     print("No image input.")

        if not self.candidate_expressions.empty:
            prompt = prompt.replace("{functions}", str(self.candidate_expressions.iloc[0:5]))

        if self.num_exprs == 1:
            example_str = r'{"eq": "c*x_0 + c",  "dim": 1}'
        elif self.num_exprs == 2:
            example_str = r'{"eq": "c*x_0 + c*x_1 + c | c*x_0 - c*x_1",  "dim": 2}'
        elif self.num_exprs == 3:
            example_str = r'{"eq": "c*x_0 - c*x_1 + c*x_2 + c | c*x_1*x_2 - c*x_0 | c*x_0 + c*x_1 - c*x_2",  "dim": 3}'
        elif self.num_exprs == 4:
            example_str = r'{"eq": "c*x_0*x_1 + c*x_2 - c*x_3 + c | c*x_1/x_0 + c*x_3 | c*x_2 - c*x_0*x_3 | c*x_3 + c*x_0 - c*x_1",  "dim": 4}'
        # import ipdb; ipdb.set_trace()
        prompt = prompt.replace("{example_str}", example_str)

        self.current_prompt = prompt

        return prompt

        
    
    def generate(self, exprs_list, init_conds_list):
        """
        this code is adopeted from odebench
        """
        solutions = []
        for i, exprs in tqdm(enumerate(exprs_list)):
            var_symbols = self.dataset.var_symbols
            # import ipdb; ipdb.set_trace()
            callable_fn = lambda t, x: np.array([f(*x) for f in [sympy.lambdify(var_symbols, eq, 'numpy') for eq in exprs]])
            init_conds = init_conds_list[i]
            # suppress invalid-value-in-power warnings during integration
            with np.errstate(invalid='ignore', divide='ignore'):
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
        import re

        # Extract all dict-like blocks
        dict_blocks = re.findall(r'\{[^}]*\}', str_funcs)
        all_expressions = []

        for block in dict_blocks:
            # Extract the equation string after "eq":
            eq_match = re.search(r'"eq"\s*:\s*"([^"]+)"', block)
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

        return all_expressions
    

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


    def update_candidates(self, expr_list, score=None):
        """
        Update the candidate expressions with the new expression.
        """
        # import ipdb; ipdb.set_trace()
        for expr in expr_list:
            if not self.check_duplicate(expr):
                # import ipdb; ipdb.set_trace()
                fit_score = self.fit_score(expr)
                if fit_score is None:
                    self.logger.info(f"Expression {expr} is not valid.")
                    continue
                # import ipdb; ipdb.set_trace()
                tmp  ={"expr": str(expr)}
                tmp.update(**fit_score)
                self.candidate_expressions = pd.concat([self.candidate_expressions, pd.DataFrame({k: [v] for k, v in tmp.items()})], ignore_index=True)

                # import ipdb; ipdb.set_trace()
            else:
                self.logger.info(f"Expression {expr} already exists in the candidate expressions.")
        # import ipdb; ipdb.set_trace()


        self.candidate_expressions = (
            self.candidate_expressions
            .drop_duplicates(subset=['expr'])
            .sort_values(by="R2", ascending=False)
            .reset_index(drop=True)
            .iloc[:100]
        )
        # import ipdb; ipdb.set_trace()
        return self.candidate_expressions
    
    def fit_score(self, exprs):
        """
        Score the expressions based on the number of variables and equations.
        """
        # import ipdb; ipdb.set_trace()
        try:
            fit_exprs = self.fit(exprs, self.init_conds, self.time_series)
        except Exception as e:
            print(f"Error in fitting: {e}")
            return None
        # import ipdb; ipdb.set_trace()
        try:
            # import ipdb; ipdb.set_trace()
            pred_data = self.generate([fit_exprs], [self.init_conds])[0]
            r2 = utils.r2_score(self.time_series, pred_data)
            complexity = [e.count_ops() for e in exprs]
        except Exception as e:
            print(f"Error in generating data: {e}")
            return None
    
        return {"fit_exprs": str(fit_exprs) , "R2": r2, "Complexity": np.sum(complexity).item()}

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


    def fit(self, raw_exprs, init_conds, gt_data):
        import numpy as np, sympy
        from scipy.integrate import solve_ivp
        from scipy.optimize import minimize
        # replace constants with symbols
        exprs_c, c_syms = self.replace_coefficients(raw_exprs)
        var_syms = self.dataset.var_symbols
        # import ipdb; ipdb.set_trace()
        # lambdify with numpy backend
        fns = [sympy.lambdify(var_syms + c_syms, expr, "numpy") for expr in exprs_c]
        # prepare time and data
        t = np.array(self.solve_config["t_eval"])
        X = gt_data if gt_data.shape[1] == len(var_syms) else gt_data[:,1:]
        # objective: negative variance-weighted R²
        def obj(c_vals):
            def rhs(tt, yy):
                return np.array([fn(*yy, *c_vals) for fn in fns], dtype=float)
            with np.errstate(invalid='ignore', divide='ignore'):
                sol = solve_ivp(rhs, (t[0], t[-1]), init_conds, t_eval=t, **{k:v for k,v in self.solve_config.items() if k not in ("t_span","t_eval")})
            if not sol.success or sol.y.shape[1] != len(t):
                return 1e6
            Y = sol.y.T
            return -utils.r2_score(X, Y)
        # no constants case
        if not c_syms:
            return exprs_c
        # initial guess zeros
        x0 = np.zeros(len(c_syms))

        # import ipdb; ipdb.set_trace()

        max_retry = 2
        for i in range(max_retry):
            # randomize initial guess
            if i > 0:
                x0 = np.random.randn(len(c_syms))*2
            # run BFGS
            res = minimize(obj, x0, method="BFGS", options={"gtol":1e-4, "maxiter":200})
            if res.success:
                break
        import ipdb; ipdb.set_trace()             
        if not res.success:
            self.logger.warning(f"fit_constants did not converge: {res.message}")
        # substitute optimized constants back
        # import ipdb; ipdb.set_trace()
        fitted = []
        for expr in exprs_c:
            expr = expr.subs({c_syms[i]: float(res.x[i]) for i in range(len(c_syms))})
            fitted.append(expr)
        return fitted
    


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
        return eq

