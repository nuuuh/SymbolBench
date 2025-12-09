import os
import re
import numpy as np
import pandas as pd

from torch import dtype, device
from torch.cuda import get_device_name, is_available
from collections.abc import Callable
from typing import Tuple, List, Dict, Any

from sklearn.metrics import mean_squared_error, r2_score
# from models import HuggingFaceModel, OpenAIModel

import sympy
from sympy.parsing.sympy_parser import parse_expr




def compute_ned(pred_exprs, true_exprs, var_symbols):
    # compute tree edit distances (Zhang-Shasha) between expressions
    from zss import simple_distance, Node
    import sympy as sp
    # helper: drop constant part

    pred_exprs = [sp.sympify(expr) for expr in pred_exprs]
    true_exprs = [sp.sympify(expr) for expr in true_exprs]

    def drop_const(expr):
        _, core = expr.as_independent(*var_symbols)
        return core
    # convert expression to tree
    def to_tree(expr):
        e = drop_const(expr)
        if not e.args:
            return Node(str(e))
        node = Node(type(e).__name__)
        for arg in e.args:
            node.addkid(to_tree(arg))
        return node
    def count_nodes(node):
        return 1 + sum(count_nodes(c) for c in node.children)
    # helper: pretty-print tree labels
    def print_tree(node, indent=0):
        print('  ' * indent + str(node.label))
        for child in node.children:
            print_tree(child, indent + 1)
    edit_distance= []
    for pe, te in zip(pred_exprs, true_exprs):
        t1, t2 = to_tree(pe), to_tree(te)
        # print("Predicted expression tree:")
        # print_tree(t1)
        # print("True expression tree:")
        # print_tree(t2)
        d = simple_distance(t1, t2)
        edit_distance.append(d)
    return edit_distance




def to_serializable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [to_serializable(v) for v in obj]
    else:
        return obj


def get_job_id() -> str:
    """
    Gets the SLURM job id from the environment variables if available.
    
    Returns
    -------
    job_id -> The SLURM job id (or None if not available).
    """
    
    job_id = os.environ.get("SLURM_JOB_ID", None)
    if job_id is not None:
        job_id += "_" + os.environ.get("SLURM_ARRAY_TASK_ID", None) if "SLURM_ARRAY_TASK_ID" in os.environ else ""
        
    return job_id

# def load_model(model_name: str, device: device, dtype: dtype, cache_dir: str = None, model_args = None) -> Any:
#     """
#     Utility to load a model from the HuggingFace model hub.
#     Mostly needed to deal with LLaVA models, that are not available on the model hub yet.

#     Parameters
#     ----------
#     model_name -> the name of the model to load.
#     device -> the device to load the model on.
#     dtype -> the dtype to load the model with.
#     cache_dir -> the cache directory to use for the model.

#     Returns
#     -------
#     model -> the loaded model.
#     """ 
#     if 'llava' in model_name:
#         model = LLaVaModelHF(model_name, device, dtype, cache_dir, **model_args)
#     elif 'gpt' in model_name:
#         model = OpenAIModel(model_name, device, dtype, cache_dir, **model_args)
#     else:
#         model = HuggingFaceModel(model_name, device, dtype, cache_dir, **model_args)

#     return model

def get_messages(prompt: str, splits: List[str] = ["system", "user"]) -> List[Dict[str, str]]:
    """
    Converts a prompt string into a list of messages for each split.
    
    Parameters:
        prompt (str): The prompt string.
        splits (list[str]): A list of the splits to parse. Defaults to ["system", "user"].
        
    Returns:
        list[dict[str, str]]: A dictionary of the messages for each split.
    """
    
    messages = []
    for split in splits:
        start_tag = f"<{split}>"
        end_tag = f"</{split}>"

        start_idx = prompt.find(start_tag)
        end_idx = prompt.find(end_tag)
        
        # Skip if the split is not in the prompt (e.g. no system prompt)
        if start_idx == -1 or end_idx == -1:
            continue
        messages.append({
            "role": split,
            "content": prompt[start_idx + len(start_tag):end_idx].strip()
        })
    
    # If no splits at all, assume the whole prompt is a user message
    if len(messages) == 0:
        messages.append({
            "role": "user",
            "content": prompt
        })

    return messages

def load_points(file_path: str) -> np.ndarray:
    """
    Loads a set of points from a file.

    Parameters
    ----------
    file_path -> the path to the file containing the points.

    Returns
    -------
    points -> the points.
    """
    if file_path.endswith(".npy"):
        points = np.load(file_path)
    elif file_path.endswith(".txt"):
        points = np.loadtxt(file_path)
    elif file_path.endswith(".csv"):
        points = pd.read_csv(file_path).values
    elif file_path.endswith(".tsv"):
        points = pd.read_csv(file_path, sep="\t").values
    else:
        raise ValueError("Invalid file format. (only .npy, .txt, .csv, and .tsv are supported)")
    return points

def normalize_points(points: np.ndarray, method: str = "minmax", percentile: int = None) -> np.ndarray:
    """
    Normalizes a set of points.

    Parameters
    ----------
    points -> the points to normalize.
    method -> the normalization method to use. (minmax, zscore, percentile)
    percentile -> the percentile to use for percentile normalization (if applicable).

    Returns
    -------
    points -> the normalized points.
    """
    if method == "percentile" and percentile is None:
        raise ValueError("Percentile normalization requires a percentile value.")
    
    ys = np.array([point[-1] for point in points])
    if method == "minmax":
        points = np.array([np.concatenate([point[:-1], [(y - ys.min()) / (ys.max() - ys.min())]]) for point, y in zip(points, ys)])
    elif method == "zscore":
        points = np.array([np.concatenate([point[:-1], [(y - ys.mean()) / ys.std()]]) for point, y in zip(points, ys)])
    elif method == "percentile":
        points = np.array([np.concatenate([point[:-1], [y /np.percentile(ys, percentile)]]) for point, y in zip(points, ys)])
    else:
        raise ValueError("Invalid normalization method.")

    points = np.round(points, 4)
    return points

def decimate_points(points: np.ndarray, max_points: int) -> np.ndarray:
    """
    Reduces the number of points to a maximum number to be used in the prompt.
    
    Parameters
    ----------
    points -> the points to decimate.
    max_points -> the maximum number of points to keep.
    
    Returns
    -------
    points -> the decimated points.
    """
    
    if points.shape[0] <= max_points:
        return points
    
    # Find an evenly spaced subset of points
    indices = np.linspace(0, points.shape[0] - 1, max_points, dtype=int)
    points = points[indices]
    return points
    
def split_points(points: np.ndarray, test_fraction: float, split_strategy: str = "random", seed: int = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Splits a set of points into train and test sets.

    Parameters
    ----------
    points -> the points to split.
    test_fraction -> the fraction of points to use for the test set.
    split_strategy -> the strategy to use for splitting the points. (random, middle, end)
    seed -> the seed to use for the random split.

    Returns
    -------
    train_points -> the train points.
    test_points -> the test points.
    """
    num_points = points.shape[0]
    num_test_points = int(num_points * test_fraction)
    points = points[points[:, 0].argsort()]

    if seed is not None:
        np.random.seed(seed)

    #! Middle and end are not working properly with n_variables > 1, not fixed as unused in final version
    if split_strategy == "random":
        indices = np.random.choice(num_points, num_test_points, replace=False)
        mask = np.ones(num_points, dtype=bool)
        mask[indices] = False
        train_points = points[mask]
        test_points = points[~mask]
    elif split_strategy == "middle":
        start = (num_points - num_test_points) // 2
        end = start + num_test_points
        train_points = np.concatenate([points[:start], points[end:]])
        test_points = points[start:end]
    elif split_strategy == "end":
        train_points = points[:-num_test_points]
        test_points = points[-num_test_points:]
    else:
        raise ValueError("Invalid split strategy.")

    return train_points, test_points

def array_to_string(points: np.ndarray, func_type='standard') -> str:
    """
    Converts a numpy array of points to a string.

    Parameters
    ----------
    points -> the numpy array of points to convert.

    Returns
    -------
    points -> the string of points.
    """
    if func_type == "DE":
        points = points[:, 1:]
    points = points.tolist()
    points_str = ""
    for point in points:
        point_str = ", ".join([str(np.round(x, 2)) for x in point])
        point_str = f"({point_str})"
        points_str += point_str + ", "

    return points_str[:-2]

def string_to_array(points: str) -> np.ndarray:
    """
    Converts a string of points to a numpy array.

    Parameters
    ----------
    points -> the string of points to convert.

    Returns
    -------
    points -> the numpy array of points.
    """
    points = points.replace("(", "").split("), ")
    points = [point.replace(")", "") for point in points]
    points = [point.split(", ") for point in points]
    points = [[float(coordinate) for coordinate in point] for point in points]
    return np.array(points)

def eval_function(function: sympy.core.function.Function, Xs: np.ndarray, num_variables: int) -> float:
    """
    Evaluates a sympy function at a point.

    Parameters
    ----------
    function -> the function to evaluate.
    Xs -> the points to evaluate the function at. (Variables have to be sorted alphabetically)
    num_variables -> the number of variables the function takes.
    
    Returns
    -------
    ys -> the value of the function at x.
    """
    symbols = function.free_symbols
    symbols = sorted(symbols, key=lambda x: str(x))
    if Xs.shape[-1] != num_variables:
        Xs = np.array(list(zip(*[x.flat for x in Xs])))
    # import ipdb; ipdb.set_trace()
    ys = []
    for point in Xs:
        if type(point) == np.ndarray:
            subs = {symbol: value for symbol, value in zip(symbols, point)}
        else:
            subs = {symbols[0]: point}
        try:
            y = function.evalf(subs=subs)
            y = float(y)
        except Exception as e:
            print(f"Error evaluating function: {function} at point {point}. {e}")
            y = np.inf
        ys.append(y)

    ys = np.array(ys)
    ys = ys.astype(np.float32)
    return ys


def eval_differential_equation(equations: List[sympy.Eq], init_conditions: Dict[str, float], 
                             x_range: np.ndarray) -> np.ndarray:
    """
    Evaluates a system of coupled differential equations using scipy's solve_ivp.

    Parameters
    ----------
    equations -> List of sympy Equations representing the coupled differential equations
                Example: For the system:
                dx/dt = y
                dy/dt = -x
                Input as [Eq(diff(x(t), t), y(t)), Eq(diff(y(t), t), -x(t))]
    
    init_conditions -> dictionary of initial conditions
                      Keys are tuples of (variable_name, derivative_order)
                      Values are the initial values
                      Example: {
                          ('x', 0): 1.0,  # x(0) = 1
                          ('y', 0): 0.0,  # y(0) = 0
                      }
    
    x_range -> numpy array of points to evaluate the solution at
               Must be sorted and include the initial condition point

    Returns
    -------
    ys -> numpy array of solution values at the given points,
          with shape (len(x_range), num_equations)
    """
    try:
        from scipy.integrate import solve_ivp
        t = sympy.Symbol('t')
        # Extract variables and their functions
        var_funcs = {}
        for eq in equations:
            for func in eq.atoms(sympy.Function):
                if str(func.func).startswith('x') and len(str(func.func)) > 1:
                    var_funcs[str(func.func)] = func
                # var_funcs[str(func.func)] = func
        
        # Create the system of first-order ODEs
        system_order = len(equations)
        y0 = []
        for i, var_name in enumerate(sorted(var_funcs.keys())):
            for order in range(1):  # For now, handling only first-order systems
                # import ipdb; ipdb.set_trace()
                try:
                    y0.append(init_conditions.get(var_name))
                except:
                    y0.append(init_conditions[i])
        # import ipdb; ipdb.set_trace()
        # Convert equations to numerical functions
        def system(t_value, y_vec):
            # Create substitution dictionary
            subs = {t: float(t_value)}
            for i, var_name in enumerate(sorted(var_funcs.keys())):
                subs[var_funcs[var_name]] = float(y_vec[i])
            # Evaluate each equation
            dydt = []
            for eq in equations:
                if isinstance(eq.lhs, sympy.Derivative):
                    rhs = eq.rhs
                else:
                    rhs = sympy.solve(eq, eq.lhs)[0]
                dydt.append(float(rhs.subs(subs)))

            return dydt
        # import ipdb; ipdb.set_trace()
        # Solve using scipy.integrate.solve_ivp
        x_range = x_range.squeeze()
        solution = solve_ivp(
            system,
            t_span=(x_range.min(), x_range.max()),
            y0=y0,
            t_eval=x_range,
            method='RK45',
            rtol=1e-3,    # Default is 1e-3
            atol=1e-6     # Default is 1e-6
        )
        
        # Extract solution values
        ys = solution.y.T  # Transpose to get (time_points, variables)
        
    except Exception as e:
        print(f"Error solving differential equations: {equations}. {e}")
        ys = np.full((len(x_range), len(equations)), np.inf)

    return ys.astype(np.float32)


def balance_brackets(function: str) -> str:
    """Helper function to balance brackets in a mathematical expression"""
    stack = []
    result = []
    
    # First pass: remove extra closing brackets
    for char in function:
        if char == '(':
            stack.append(char)
            result.append(char)
        elif char == ')':
            if stack:  # If there's a matching opening bracket
                stack.pop()
                result.append(char)
            # Skip extra closing brackets
        else:
            result.append(char)
    
    # Add missing closing brackets at the end
    result.extend([')'] * len(stack))
    
    return ''.join(result)


def clean_function(function: str) -> str:
    """
    Cleans a function string to be evaluable.
    """
    ori_function = function
    function_list = []

    for function in ori_function.split("|"):
        function = function.strip(".")
        function = function.replace(" ", "")
        function = function.replace("cdot", "*")

        if "=" in function:
            function = function.split("=")[1]
        elif ":" in function:
            function = function.split(":")[1]

        # Remove characters that are not allowed in a function
        removals = ["'", '"', "\\", "\n", "\t", "\r", " ", "_"]
        for removal in removals:
            function = function.replace(removal, "")

        # Remove trailing operators
        while len(function) > 1 and function[-1] in ["+", "-", "*", "/", "**"]:
            if len(function) == 1:
                return lambda x: 0
            function = function[:-1]

        # Remove leading operators
        while len(function) > 1 and function[0] in ["+", "*", "/", "**"]:
            if len(function) == 1:
                return lambda x: 0
            function = function[1:]

        # Remove leading indicators of a function definition
        removals = ["Function", "Newfunction", "Thefunctionis", ":"]

        for removal in removals:
            if removal.lower() in function.lower():
                function = function.replace(removal, "")
                function = function.strip()
        function = re.sub(r'x\(x(\d+)\)', r'x\1', function)
        # Balance brackets by removing extra ones
        function = balance_brackets(function)

        function_list.append(function)

    # import ipdb; ipdb.set_trace()

    if len(function_list) == 1:
        return function_list[0]
    else:
        return ' | '.join(function_list)

def string_to_function(function: str, num_variables: int = 1, func_type="standard") -> Callable[[float], float]:
    """
    Converts a string to a callable function using eval.

    Parameters
    ----------
    function -> the string to convert.
    num_variables -> the number of variables the function should take.

    Returns
    -------
    f -> the callable function.
    """
    function = clean_function(function)
    if func_type == "standard":
        np_func = ["sin", "cos", "tan", "exp", "log", "sqrt"]
        function = function.replace("^", "**")
        #! This only works for variables in x (x1, x2, x3, ...)
        #! This only works with coefficients that end with numbers (e.g. c0, c1, c2, ...)
        function = re.sub(r"(\d)x", r"\1*x", function)
        regex = r"(\d)(" + "|".join(np_func) + ")"
        function = re.sub(regex, r"\1*\2", function)
        f = parse_expr(function)
        return f
    
    if func_type == "DE":
        # Split equations by |
        eqs = [eq.strip() for eq in function.split('|')]
        num_eqs = len(eqs)
        # Define symbols
        t = sympy.Symbol('t')
        variables = []
        equations = []
        
        # Create function symbols (x_0(t), x_1(t), etc)
        for i in range(1, len(eqs)+1):
            var = sympy.Function(f'x{i}')(t)
            variables.append(var)
        # Parse each equation
        for i, eq_str in enumerate(eqs):
            # Replace x_0, x_1 with function calls
            for j, var in enumerate(variables):
                # eq_str = eq_str.replace(f'x{j+1}', str(var)
                var_name = f'x{j+1}'
                eq_str = re.sub(f'{var_name}(?!\\(t\\))', str(var), eq_str)
            
            # # Replace c_0, c_1, etc with symbols
            # for match in re.finditer(r'c(\d+)', eq_str):
            #     c_idx = match.group(1)
            #     eq_str = eq_str.replace(f'c{c_idx}', f'c{c_idx}')
            
            np_func = ["sin", "cos", "tan", "exp", "log", "sqrt"]
            eq_str = eq_str.replace("^", "**")
            #! This only works for variables in x (x1, x2, x3, ...)
            #! This only works with coefficients that end with numbers (e.g. c0, c1, c2, ...)
            eq_str = re.sub(r"(\d)x", r"\1*x", eq_str)
            regex = r"(\d)(" + "|".join(np_func) + ")"
            eq_str = re.sub(regex, r"\1*\2", eq_str)
            
            # import ipdb; ipdb.set_trace()
            # Parse the right-hand side
            rhs = parse_expr(eq_str)

            
            # Create the equation
            eq = sympy.Eq(variables[i].diff(t), rhs)
            equations.append(eq)
        
        assert len(equations) == num_eqs, f"Number of equations {len(equations)} does not match expected {num_eqs}."

        return equations

    




def is_valid_function(function: str, current_functions: Any, num_variables: int = 1, func_type='standard', num_eqs=1) -> Tuple[bool, str]:
    """
    Checks if a function is valid.

    Parameters
    ----------
    function -> the function to check.
    current_functions -> the current functions in the prompt.
    num_variables -> the number of variables the function should take.

    Returns
    -------
    valid -> whether the function is valid.
    reason -> the reason the function is invalid (if applicable).
    """ 
    valid = True
    reason = ""

    # Check for invalid x() pattern in the string representation
    if isinstance(function, str):
        invalid_pattern = r'x\([^0-9].*?\)|x\(\)'
        if re.search(invalid_pattern, function):
            valid = False
            reason += "Invalid variable syntax: x() must use numerical indices like x1, x2, etc. "
            return valid, reason

    if type(function) == str:
        # import ipdb; ipdb.set_trace()
        f_list = string_to_function(function, num_variables, func_type)
        if type(f_list) != list:
            f_list = [f_list]
    else:
        f_list = [function]
    
    if len(f_list) != num_eqs:
        valid = False
        reason += f"Invalid number of equations: expected {num_eqs}, got {len(f_list)}. "
        return valid, reason

    # import ipdb; ipdb.set_trace()
    for f in f_list:
        symbols = f.free_symbols
        variables = [str(symbol) for symbol in symbols if str(symbol).startswith("x") or str(symbol).startswith("t")]

        if len(variables) > num_variables:
            valid = False
            reason += "Too many variables in function."

        if current_functions is not None and current_functions.func_in_list(f):
            valid = False
            reason += "Function already in prompt."

    return valid, reason

def format_exp(x: float, d: int = 6) -> str:
    """
    Formats a number in scientific notation with custom precision. (used in Scorers)

    Parameters
    ----------
    x -> the number to format.
    d -> the number of decimal places to round to.

    Returns
    -------
    x -> the formatted number.
    """
    n = int(np.floor(np.log10(abs(x))))
    significand = x / 10 ** n
    exp_sign = '+' if n >= 0 else '-'
    return f'{significand:.{d}f}e{exp_sign}{n:02d}'

# def func_equals(f1: Any, f2: Any, num_variables: int) -> bool:
#     """
#     Checks if two functions are equal. Used in place of sympy.equals as the latter can become very slow for certain functions.
#     https://stackoverflow.com/questions/37112738/sympy-comparing-expressions

#     Parameters
#     ----------
#     f1 -> the first function.
#     f2 -> the second function.
#     num_variables -> the number of variables the functions should take.

#     Returns
#     -------
#     equal -> whether the functions are equal.
#     """
#     if f1 == f2:
#         return True
#     if f1 is None or f2 is None:
#         return False
#     if f1.free_symbols != f2.free_symbols:
#         return False
#     if f1.free_symbols != set([sympy.Symbol(f"x{i + 1}") for i in range(num_variables)]):
#         return False
#     return False


def func_equals(f1: Any, f2: Any, num_variables: int) -> bool:
    """
    Checks if two functions or sets of functions are equal.
    Used in place of sympy.equals as the latter can become very slow for certain functions.
    Handles both single functions and coupled differential equations (lists/tuples).
    
    Parameters
    ----------
    f1 -> the first function or set of functions.
    f2 -> the second function or set of functions.
    num_variables -> the number of variables the functions should take.
    
    Returns
    -------
    equal -> whether the functions are equal.
    """
    # Handle lists or tuples of functions (coupled differential equations)
    if isinstance(f1, (list, tuple)) and isinstance(f2, (list, tuple)):
        # Must have same number of equations
        if len(f1) != len(f2):
            return False
            
        # Check for direct equality
        if f1 == f2:
            return True
            
        # Try all possible matchings for unordered comparisons
        from itertools import permutations
        for perm in permutations(f2):
            all_equal = True
            for eq1, eq2 in zip(f1, perm):
                if not func_equals(eq1, eq2, num_variables):
                    all_equal = False
                    break
            if all_equal:
                return True
        return False
        
    # From here on, we're handling single functions
    
    # Direct equality check (fast path)
    if f1 == f2:
        return True
        
    # Basic checks
    if f1 is None or f2 is None:
        return False
        
    # Check for symbolic equations (sympy.Eq objects)
    if hasattr(f1, 'lhs') and hasattr(f2, 'lhs'):
        # For differential equations, compare both sides
        lhs_equal = func_equals(f1.lhs, f2.lhs, num_variables)
        rhs_equal = func_equals(f1.rhs, f2.rhs, num_variables)
        return lhs_equal and rhs_equal
    
    # For sympy expressions, check free symbols
    if hasattr(f1, 'free_symbols') and hasattr(f2, 'free_symbols'):
        if f1.free_symbols != f2.free_symbols:
            return False
            
        # For standard functions, we expect symbols like x1, x2, etc.
        if all(str(s).startswith('x') and not str(s).startswith('x(') for s in f1.free_symbols):
            # Sample random points and compare function values
            import numpy as np
            
            # Generate test points
            n_tests = 10
            test_points = np.random.uniform(-5, 5, (n_tests, num_variables))
            
            # Convert symbols to ordered list for substitution
            symbols1 = sorted(list(f1.free_symbols), key=lambda x: str(x))
            symbols2 = sorted(list(f2.free_symbols), key=lambda x: str(x))
            
            # Test at multiple points
            for point in test_points:
                try:
                    # Create substitution dictionaries
                    subs1 = {s: float(val) for s, val in zip(symbols1, point)}
                    subs2 = {s: float(val) for s, val in zip(symbols2, point)}
                    
                    # Evaluate both functions
                    val1 = float(f1.evalf(subs=subs1))
                    val2 = float(f2.evalf(subs=subs2))
                    
                    # If values differ significantly, functions are not equal
                    if abs(val1 - val2) > 1e-6:
                        return False
                except Exception:
                    # If evaluation fails, assume they're different
                    return False
            
            # If all test points match, consider them equal
            return True
    
    # Default: not equal
    return False



def count_nodes(formula: Any) -> int:
    """
    Gets the complexity of a sympy formula, represented by the number of nodes in its expression tree.
    
    Parameters
    ----------
    formula -> the formula to get the complexity of.
    
    Returns
    -------
    complexity -> the complexity of the formula.
    """
    if type(formula) is list or type(formula) is tuple:
        complexity = 0
        for f in formula:
            complexity += f.count_ops()
        return complexity
    else:
        return formula.count_ops()

def replace_zero_coefficients(expr: Any, formula: Any, threshold: float = 1e-2) -> Any:
    """
    Replaces coefficients that are close to zero in a formula with zero.
    
    Parameters
    ----------
    expr -> the expression to replace coefficients in (with coefficients c0, c1...)
    formula -> the formula to replace coefficients in (with numerical coefficients)
    threshold -> the threshold to consider a coefficient zero.
    
    Returns
    -------
    expr -> the expression with zero coefficients replaced.
    formula -> the formula with zero coefficients replaced.
    """
    coeffs_dict = formula.as_coefficients_dict()
    expr_dict = expr.as_coefficients_dict()
    
    for key, value in coeffs_dict.items():
        if abs(value) < threshold:
            expr_dict[key] = 0
            formula = formula.subs(key, 0)
            
    expr = expr.subs(expr_dict)
    
    print(expr, formula)