import json
from odeformer.odebench.strogatz_equations import equations
from odeformer.odebench.solve_and_plot import config, process_equations, solve_equations, plot_prediction
import pandas as pd
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from odeformer.metrics import r2_score
from tqdm import tqdm 
import torch
import os
from sympy import sympify

from utils.utils import compute_ned

import argparse
import time  # for timeout mechanism
parser = argparse.ArgumentParser(description='Evaluate DE models')
parser.add_argument('--dataset', type=str, default='strogatz_extended')
parser.add_argument('--model', type=str, default='odeformer', help='Model to evaluate')
parser.add_argument('--eval_type', type=str, default='', help='Evaluation type')
parser.add_argument('--OOD', type=int, default=0, help='Evaluate on OOD data')
parser.add_argument('--timeout', type=int, default=15, help='Timeout for ODE solving in seconds')
parser.add_argument('--folder', type=str, default='')
parser.add_argument('--LLM', type=int, default=0)
parser.add_argument('--modality', type=str, default='textual')
args = parser.parse_args()
TIMEOUT = args.timeout


if args.dataset == 'strogatz_extended' or args.dataset == "strogatz":
    if os.path.exists('data/gt_DEs.pt'):
        gt_equations = torch.load('data/gt_DEs.pt', weights_only=False)
    else:
        gt_equations = equations
        process_equations(gt_equations)
        solve_equations(gt_equations, config)
        torch.save(gt_equations, 'data/gt_DEs.pt')
    dataset = args.dataset
    idx_list = list(range(1, 64))

elif args.dataset == 'physiome':
    with open("data/Physiome/idx_DE", "r") as idx_f:
        idx_list = idx_f.read().strip().split()
        idx_list = [int(i) for i in idx_list]
    # import ipdb; ipdb.set_trace()
    with open('data/Physiome/solved_small_odes.json', 'r') as f:
        gt_equations = json.load(f)
    
    dataset = "physiome"

# import ipdb; ipdb.set_trace()

models=(args.model,)

eval_type = args.eval_type
modality = args.modality

solve_config = {
    "t_span": (0, 10),  # time span for integration
    "method": "LSODA",  # method for integration
    "rtol": 1e-5,  # relative tolerance (let's be strict)
    "atol": 1e-7,  # absolute tolerance (let's be strict)
    "first_step": 1e-6,  # initial step size (let's be strict)
    "t_eval": np.linspace(0, 10, 150),  # output times for the solution
    "min_step": 1e-10,  # minimum step size (only for LSODA)
}

def solve_ivp_with_timeout(fun, init_cond):
    """Run solve_ivp with time-wrapper to enforce timeout."""
    start_time = time.time()
    def fun_wrapper(t, x):
        if time.time() - start_time > TIMEOUT:
            raise RuntimeError(f"ODE solve timed out after {TIMEOUT} seconds")
        return fun(t, x)
    sol = solve_ivp(fun_wrapper, **solve_config, y0=init_cond)
    return sol

OOD = args.OOD

for model in models:
    if len(args.folder) > 0:
        folder = args.folder
    else:
        folder = None
    if model == "odeformer":
        data_path = f"baseline_results/{dataset}/odeformer/eval_{dataset}.csv"
    elif args.LLM:
        if folder is not None:
            data_path = os.path.join(folder, f"eval_results.csv")
        else:
            data_path = f"runs/{model}/{args.modality}/{eval_type.split('_')[-1]}_{dataset}/eval_results.csv"
    else:
        data_path = f"baseline_results/{dataset}/{model}/eval_result.csv"
    pred_result = pd.read_csv(data_path)
    # pred_result['info_name'] = pred_result['info_name'].apply(lambda x: x[:-3])

    success = 0
    results = []
    for i, gt_func in tqdm(enumerate(gt_equations)):
        idx = gt_func['id']
        if idx not in idx_list:
            continue
        try:
            gt_eq_str = gt_func['substituted']
            
            # import ipdb; ipdb.set_trace()
            try:
                tmp_result = pred_result[idx == pred_result['ids']].copy()
            except:
                try:
                    tmp_result = pred_result[idx == pred_result['index']].copy()
                except:
                    tmp_result = pred_result[idx == pred_result['id']].copy()

            try:
                eq_str = tmp_result['predicted_trees'].iloc[0]
            except:
                # import ipdb; ipdb.set_trace()
                try:
                    eq_str = tmp_result['predicted_equation'].iloc[0]
                except:
                    print("No predicted equation found")
                    continue

            # import ipdb; ipdb.set_trace()
            duration = tmp_result['duration_fit'].iloc[0]
            eqs = [e.strip() for e in eq_str.split('|')]
            dim = len(eqs)
            var_symbols = sp.symbols([f'x_{i}' for i in range(dim)])
            exprs = [sp.sympify(eq) for eq in eqs]
            # ensure constants (e.g. pi, E) are numeric floats
            exprs = [expr.evalf() for expr in exprs]
            fns = [sp.lambdify(var_symbols, expr, 'numpy') for expr in exprs]

            # gt_eqs = [e.strip() for e in gt_eq_str.split('|')]
            gt_eqs = gt_eq_str[0]
            gt_exprs = [sp.sympify(eq) for eq in gt_eqs]
            gt_complexity = sum(expr.count_ops() for expr in gt_exprs)
            gt_funcs_sympy = [sp.lambdify(var_symbols, expr, 'numpy') for expr in gt_exprs]

            def ode_func(t, x):
                return np.array([f(*x) for f in fns])
            # import ipdb; ipdb.set_trace()
            def gt_ode(t, x):
                return np.array([f(*x) for f in gt_funcs_sympy])

            initial_conditions = gt_func['init']
            # 0 for ID, 1 for OOD

            if OOD:
                if len(initial_conditions) == 1:
                    id_init = initial_conditions[0]
                    init_cond = np.abs(id_init + np.random.normal(0, 0.1, size=len(id_init)))
                    # import ipdb; ipdb.set_trace()
                    sol = solve_ivp_with_timeout(gt_ode, init_cond)
                    gt_tra = sol.y.T

                else:    
                    init_cond = initial_conditions[1]
                    gt_tra = np.array(gt_func['solutions'][0][1]['y']).T
            else:
                init_cond = initial_conditions[0]
                gt_tra = np.array(gt_func['solutions'][0][0]['y']).T
            # assert gt_func['solutions'][0][i]['init']  == init_cond
            print("solving for initial condition: ", init_cond)
            # sol = solve_ivp(ode_func, config['t_span'], init_cond, t_eval=config['t_eval'], method='LSODA', rtol=1e-4, atol=1e-4)
            sol = solve_ivp_with_timeout(ode_func, init_cond)
            
            pred_tra = np.array(sol.y.T)

            # import ipdb; ipdb.set_trace()
            R2 = r2_score(gt_tra, pred_tra)
            complexity = sum(expr.count_ops() for expr in exprs)
            # --- begin symbolic proximity metrics (ignore constants) ---
            # compute symbolic proximity metrics via helpers
            ned_metrics = sum(compute_ned(exprs, gt_exprs, var_symbols))
            # --- end symbolic proximity metrics ---

            # import ipdb; ipdb.set_trace()

            print(R2)
            

            results.append({
                'index': idx,
                'num_eqs': dim,
                'predicted_equation': eq_str,
                'groundtruth_equation': gt_eq_str,
                'R2': R2,
                'predicted_complexity': complexity,
                'groundtruth_complexity': gt_complexity,
                'initial_conditions': init_cond,
                "duration_fit": duration,
                'ned': ned_metrics,
                # 'groundtruth_causal_graph': gt_causal_graph,
            })
            success+= 1
        except Exception as e:
            results.append({
                    'index': idx,
                    'num_eqs': None,
                    'predicted_equation': None,
                    'groundtruth_equation': None,
                    'R2': None,
                    'predicted_complexity': None,
                    'groundtruth_complexity': None,
                    'initial_conditions': None,
                    "duration_fit": None,
                    'ned': None,
                })
            print(f"Error processing equation at index {idx}: {e}")

    print(success, "successful equations out of ", len(gt_equations))
    # Save results to DataFrame
    results_df = pd.DataFrame(results)
    if OOD:
        if len(eval_type) > 0:
            results_df.to_csv(f"baseline_results/DE_eval_OOD/{model}_{modality}_{dataset}_{eval_type}.csv", index=False)
        else:
            results_df.to_csv(f"baseline_results/DE_eval_OOD/{model}_{modality}_{dataset}.csv", index=False)
    else:
        if len(eval_type) > 0:
            results_df.to_csv(f"baseline_results/DE_eval_ID/{model}_{modality}_{dataset}_{eval_type}.csv", index=False)
        else:
            results_df.to_csv(f"baseline_results/DE_eval_ID/{model}_{modality}_{dataset}.csv", index=False)

