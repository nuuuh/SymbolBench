from odeformer.odebench.strogatz_equations import equations
from odeformer.odebench.solve_and_plot import config, process_equations, solve_equations, plot_prediction
import pandas as pd
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from odeformer.metrics import r2_score
from tqdm import tqdm 
import networkx as nx
import matplotlib.pyplot as plt
import io
import base64
import os
import torch
import json


def extract_causal_graph(gt_exprs):
    """Extracts a causal graph from a list of sympy expressions."""
    variables = [str(s) for s in sorted({s for expr in gt_exprs for s in expr.free_symbols}, key=lambda x: str(x))]
    print("variables: ", variables)
    graph = {v: set() for v in variables}
    for lhs_idx, expr in enumerate(gt_exprs):
        lhs_var = variables[lhs_idx] if lhs_idx < len(variables) else None
        rhs_vars = {str(s) for s in expr.free_symbols}
        if lhs_var:
            rhs_vars.discard(lhs_var)
            graph[lhs_var] = rhs_vars
    return graph


def extract_causal_graph_with_lag(gt_exprs, max_lag=1):
    """Extracts a causal graph with lag consideration from sympy expressions."""
    variables = [str(s) for s in sorted({s for expr in gt_exprs for s in expr.free_symbols}, key=lambda x: str(x))]
    graph = {v: set() for v in variables}
    for lhs_idx, expr in enumerate(gt_exprs):
        lhs_var = variables[lhs_idx] if lhs_idx < len(variables) else None
        rhs_vars = {str(s) for s in expr.free_symbols}
        if lhs_var:
            rhs_vars.discard(lhs_var)
            # Add lagged edges
            for rv in rhs_vars:
                for lag in range(0, max_lag+1):
                    if lag == 0:
                        graph[lhs_var].add(rv)
                    else:
                        graph[lhs_var].add(f"{rv}(t-{lag})")
    return graph


def extract_causal_graph_with_lag_only_vars(gt_exprs, var_prefix='x_', max_lag=1):
    """Extracts a causal graph with lag consideration, only for variables like x_0, x_1, ..."""
    variables = [str(s) for s in sorted({s for expr in gt_exprs for s in expr.free_symbols if str(s).startswith(var_prefix)}, key=lambda x: str(x))]
    graph = {v: set() for v in variables}
    for lhs_idx, expr in enumerate(gt_exprs):
        lhs_var = variables[lhs_idx] if lhs_idx < len(variables) else None
        rhs_vars = {str(s) for s in expr.free_symbols if str(s).startswith(var_prefix)}
        if lhs_var:
            rhs_vars.discard(lhs_var)
            for rv in rhs_vars:
                for lag in range(0, max_lag+1):
                    if lag == 0:
                        graph[lhs_var].add(rv)
                    else:
                        graph[lhs_var].add(f"{rv}(t-{lag})")
    return graph


def plot_causal_graph(graph, title=None):
    G = nx.DiGraph()
    for src, tgts in graph.items():
        for tgt in tgts:
            G.add_edge(tgt, src)
    plt.figure(figsize=(4, 3))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=1200, font_size=10, arrowsize=20)
    if title:
        plt.title(title)
    buf = io.BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format='png')
    plt.close()
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode('utf-8')
    return img_b64


# if os.path.exists('data/gt_DEs.pt'):
#     gt_equations = torch.load('data/gt_DEs.pt')
# else:
#     gt_equations = equations
#     process_equations(gt_equations)
#     solve_equations(gt_equations, config)

gt_equations = equations
json_path = 'data/Physiome/aligned_small_odes.json'
with open(json_path, 'r') as f:
    gt_equations = json.load(f)

# import ipdb; ipdb.set_trace()
gt_equations = gt_equations['equations']

process_equations(gt_equations)
# import ipdb; ipdb.set_trace()
 
# torch.save(gt_equations, 'gt_DEs_inter.pt')

config = {
    "t_span": (0, 10),  # time span for integration
    "method": "LSODA",  # method for integration
    "rtol": 1e-5,  # relative tolerance (let's be strict)
    "atol": 1e-7,  # absolute tolerance (let's be strict)
    "first_step": 1e-6,  # initial step size (let's be strict)
    "t_eval": np.linspace(0, 10, 150),  # output times for the solution
    "min_step": 1e-10,  # minimum step size (only for LSODA)
}

solve_equations(gt_equations, config)

# import ipdb; ipdb.set_trace()

def remove_unpickleable(obj):
    if isinstance(obj, dict):
        return {k: remove_unpickleable(v) for k, v in obj.items() if k not in ['lambdified', 'lambdified_func']}
    elif isinstance(obj, list):
        return [remove_unpickleable(v) for v in obj]
    else:
        return obj

for eq in gt_equations:
    eq['substituted'][0] = [ str(t) for t in eq['substituted'][0]]
import ipdb; ipdb.set_trace()

gt_equations_clean = remove_unpickleable(gt_equations)
print(len(gt_equations_clean))
gt_equations_clean = [eq for eq in gt_equations_clean if len(eq['solutions'][0]) > 0 and eq['id']!=441]

print(len(gt_equations_clean))

with open('data/Physiome/solved_small_odes.json', 'w') as f:
    json.dump(gt_equations_clean, f)

# torch.save(gt_equations_clean, 'data/gt_DEs_new.pt')
# torch.save(gt_equations, 'data/gt_DEs.pt')

# fig_dir = 'causal_graph_figs'
# os.makedirs(fig_dir, exist_ok=True)

# results = []
# for idx, gt_func in tqdm(enumerate(gt_equations)):
#     gt_eq_str = gt_func['eq']
#     gt_eqs = [e.strip() for e in gt_eq_str.split('|')]
#     gt_exprs = [sp.sympify(eq) for eq in gt_eqs]
#     gt_complexity = sum(expr.count_ops() for expr in gt_exprs)
#     initial_conditions = gt_func['init']
#     for i, init_cond in enumerate(initial_conditions) :
#         gt_tra = np.array(gt_func['solutions'][0][i]['y']).T
#         # Only include variables like x_0, x_1, ... in the causal graph
#         causal_graph = extract_causal_graph_with_lag_only_vars(gt_exprs, var_prefix='x_', max_lag=1)
#         fig_b64 = plot_causal_graph(causal_graph, title=f'Causal Graph {idx}_{i}')
#         fig_path = os.path.join(fig_dir, f'causal_graph_{idx}_{i}.png')
#         with open(fig_path, 'wb') as f:
#             f.write(base64.b64decode(fig_b64))
#         results.append({
#             'differential_equations': gt_eqs,
#             'gt_tra': gt_tra,
#             'causal_graph': causal_graph,
#             'causal_graph_fig_path': fig_path
#         })



