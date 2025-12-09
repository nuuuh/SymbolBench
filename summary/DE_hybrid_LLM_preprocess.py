import numpy as np
import pandas as pd
import os
import sys
import json
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from copy import deepcopy
import sympy
import sympy as sp
from scipy.integrate import solve_ivp
from sklearn.metrics import r2_score
from utils.utils import compute_ned
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--eval_type', type=str, default="context_strogatz", help='Evaluation type') # hybrid_LLM_DE_w_llm context_strogatz
parser.add_argument('--path', type=str, default="runs/gpt-4o-mini/textual/") # hybrid_output/gpt-4o-mini/textual/ runs/gpt-4o-mini/textual/
args = parser.parse_args()

path = args.path
eval_type = args.eval_type
solve_config = {
            "t_span": (0, 10),
            "method": "LSODA",       # switch to BDF to avoid LSODA prints
            "rtol": 1e-6,
            "atol": 1e-6,
            "first_step": 1e-6,
            "t_eval": np.linspace(0, 10, 150),
            "dense_output": False,
        }

def final_evaluate(inferred, gt_info):
        
        gt_func = gt_info['substituted'][0]
        gt_func = [sympy.parse_expr(item) for item in gt_func]
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
            ID_initial_conditions = gt_info.init[0]
            OOD_initial_conditions = gt_info.init[1]

            def generate(init_conds):
                gt_sol = solve_ivp(gt_ode, **solve_config, y0=init_conds)
                inferred_sol = solve_ivp(inferred_ode, **solve_config, y0=init_conds)
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


if __name__ == "__main__":
    gt_data_path = "data/strogatz_extended.json"
    with open(gt_data_path, 'r') as f:
        gt_data = json.load(f)
    gt_data = pd.DataFrame(gt_data)

    sample_folder = os.path.join(path, eval_type)
    for folder in os.listdir(sample_folder):
        idx=int(folder.split('_')[-1])
        # import ipdb; ipdb.set_trace()
        if idx not in gt_data['id'].values:
            continue
        print(f"Processing folder: {folder}")
        try:
            data = np.load(os.path.join(sample_folder, folder, 'final_results.npy'), allow_pickle=True).item()
        except:
            print(f"Skipping {folder} due to missing final_results.npy")
            continue

        best = data['candidates'].iloc[0]['fit_exprs']
        gt_info = gt_data[gt_data['id'] == idx].iloc[0]

        metrics = final_evaluate(best, gt_info)
        results = {"exprs": best, "metrics": metrics}
        np.save(os.path.join(sample_folder, folder, 'final_summary.npy'), results)

    