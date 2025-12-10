import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sympy as sp

from typing import List, Tuple, Dict
from tigramite import data_processing as pp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

with open("data/Physiome/solved_odes.json", "r") as f:
    data = json.load(f)

with open("data/Physiome/generic_odes.json", "r") as f:
    context_data = json.load(f)

with open("data/strogatz_extended.json", "r") as f:
    stro_data = json.load(f)
stro_data = pd.DataFrame(stro_data)
context_data = pd.DataFrame(context_data)

# with open("data/strogatz_extended.json", "r") as f:
#     data = json.load(f)

data = pd.DataFrame(data)
data = data[data['dim'] > 2]




def generate_synthetic_var(
    N=5, T=500, L_max=2, graph_type="erdos_renyi",
    edge_prob=0.3, coef_range=(0.5, 1.0), noise_std=0.1, seed=1
) -> Tuple[np.ndarray, List[Tuple[str, str, int]]]:
    """
    Generate VAR synthetic data and ground-truth edges.
    Returns:
      data: shape (T, N)
      true_edges: list of (src, tgt, lag)
    """
    rng = np.random.RandomState(seed)
    # build skeleton
    if graph_type == "erdos_renyi":
        G = nx.erdos_renyi_graph(N, p=edge_prob, seed=seed, directed=True)
    elif graph_type == "small_world":
        k = min(N-1, max(2, int(edge_prob*(N-1)) + (int(edge_prob*(N-1))%2)))
        U = nx.watts_strogatz_graph(N, k=k, p=edge_prob, seed=seed)
        G = nx.DiGraph()
        for u, v in U.edges():
            if rng.rand() < 0.5: G.add_edge(u, v)
            else:                G.add_edge(v, u)
    elif graph_type == "scale_free":
        m = max(1, int(edge_prob*(N-1)))
        G = nx.barabasi_albert_graph(N, m=m, seed=seed).to_directed()
    else:
        raise ValueError("Unsupported graph_type")
    # weights per lag
    W = {tau: np.zeros((N, N)) for tau in range(1, L_max+1)}
    for tau in W:
        for i, j in G.edges():
            if rng.rand() < edge_prob:
                W[tau][j, i] = rng.uniform(*coef_range)
    # simulate VAR
    data = np.zeros((T, N))
    data[:L_max] = rng.randn(L_max, N)
    noise = rng.randn(T, N) * noise_std
    for t in range(L_max, T):
        for tau in range(1, L_max+1):
            data[t] += data[t-tau] @ W[tau].T
        data[t] += noise[t]
    # extract true edges
    var_names = [f"X{i}" for i in range(N)]
    true_edges: List[Tuple[str, str, int]] = []
    for tau, M in W.items():
        for i in range(N):
            for j in range(N):
                if M[j, i] != 0.0:
                    true_edges.append((var_names[i], var_names[j], tau))
    return data, true_edges


def wrap_for_tigramite(data: np.ndarray) -> pp.DataFrame:
    """Wrap numpy array to Tigramite DataFrame"""
    N = data.shape[1]
    var_names = [f"X{i}" for i in range(N)]
    return pp.DataFrame(data, var_names=var_names)


def extract_causal_graph_static_auto(
    eq_str: str,
    lag=None
) -> Tuple[List[Tuple[str, str]], Dict[str, List[str]]]:
    """
    Parse coupled equations RHS only ('|' separated) and extract
    static causal edges (no self-loops).
    """
    exprs = [e.strip() for e in eq_str.split('|')]
    n = len(exprs)
    var_names = [f"x_{i}" for i in range(n)]
    syms = sp.symbols(var_names)
    symmap = dict(zip(var_names, syms))
    edges: List[Tuple[str, str]] = []
    graph_dict: Dict[str, List[str]] = {v: [] for v in var_names}
    for tgt, expr in zip(var_names, exprs):
        edges.append((tgt, tgt, 1))  # self-loop

        tree = sp.sympify(expr, locals=symmap)
        for dep in tree.free_symbols:
            src = str(dep)
            if src in var_names and src != tgt:
                if lag is not None:
                    edges.append((src, tgt, lag))
                else:
                    edges.append((src, tgt))
                graph_dict[tgt].append(src)
    return edges, graph_dict


def evaluate_pcmci(true_edges, results, var_names, alpha=0.05):
    """
    Compare true_edges (list of (src, tgt, lag)) to PCMCI results
    and compute precision, recall, and F1.
    """
    p_mat = results['p_matrix']
    N, _, tau_max1 = p_mat.shape

    pred = set()
    for i in range(N):
        for j in range(N):
            for τ in range(1, tau_max1):
                if p_mat[j, i, τ] < alpha:
                    var_names[i] = var_names[i].replace("X", "x_")
                    var_names[j] = var_names[j].replace("X", "x_")
                    pred.add((var_names[i], var_names[j], τ))

    true_set = set(true_edges)
    tp = len(pred & true_set)
    fp = len(pred - true_set)
    fn = len(true_set - pred)

    precision = tp / (tp + fp) if tp + fp else 0.0
    recall    = tp / (tp + fn) if tp + fn else 0.0
    f1        = (2 * precision * recall / (precision + recall)
                 if precision + recall else 0.0)

    return {'precision': precision, 'recall': recall, 'f1': f1, 
            'tp': tp, 'fp': fp, 'fn': fn, 'pred_count': len(pred)}


def plot_causal_graph(
    edges: List[Tuple[str, str, int]],
    var_names: List[str],
    title: str,
    color: str,
    ax: plt.Axes
):
    """Draw directed-lag causal graph on given Axes"""
    G = nx.DiGraph()
    for src, tgt, lag in edges:
        G.add_edge(src, tgt, lag=lag)
    pos = nx.spring_layout(G, seed=42)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=color, node_size=600)
    nx.draw_networkx_labels(G, pos, ax=ax, font_weight='bold')
    nx.draw_networkx_edges(G, pos, ax=ax, arrowstyle='-|>', arrowsize=15)
    edge_labels = {(u, v): f"lag {d['lag']}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=edge_labels, font_color='gray')
    ax.set_title(title)
    ax.axis('off')


# ---- Main loop over a data DataFrame ----
if __name__ == "__main__":
    import pprint
    # assume 'data' DataFrame is defined elsewhere
    results = []

    for ids, (idx, row) in enumerate(data.iterrows()):
        eq_str = row['eq']

        true_edges, graph_dict = extract_causal_graph_static_auto(eq_str, lag=1)
        print("Edges:")
        for src, tgt, lag in true_edges:
            print(f"  {src} → {tgt} lag={lag}")
        print("\nDict form:")
        pprint.pprint(graph_dict)

        # extract and wrap time series
        time_series = np.array(row['solutions'][0][0]['y']).T
        try:
            dataframe = wrap_for_tigramite(time_series)
        except:
            print(f"Error wrapping time series: {idx}")
            continue
        var_names = dataframe.var_names
        
        # import ipdb; ipdb.set_trace()
        cur_context = context_data[context_data['id'] == idx].iloc[0]


        results.append({
            'idx': ids,
            'eq_index': idx,
            'eq_str': eq_str,
            'true_edges': true_edges,
            'graph_dict': graph_dict,
            'data': time_series,
            'variable_names': cur_context['x_mapping'],
            'domain': cur_context['domain'],
        })

    # import ipdb; ipdb.set_trace()

    np.save("data/science_sgc_physiome.npy", results)

    # import ipdb; ipdb.set_trace()
    sci_eqs = np.load("data/science_sgc.npy", allow_pickle=True)
    sci_results = []
    for i, item in enumerate(sci_eqs):
        idx = i+len(results)
        eq_index = item['index']
        stro_item = stro_data[stro_data.id == eq_index].iloc[0]
        # import ipdb; ipdb.set_trace()
        variable_names = stro_item['var_description']
        
        sci_results.append({
            'idx': idx,
            'eq_index': eq_index,
            'eq_str': item['eq_str'],
            'true_edges': item['true_edges'],
            'graph_dict': item['graph_dict'],
            'data': item['data'],
            'variable_names': variable_names,
            'domain': None,
        })

    import ipdb; ipdb.set_trace()
    combined = results + sci_results
    np.save("data/science_sgc_total.npy", combined)
        
        # run PCMCI
        # parcorr = ParCorr(significance='analytic')
        # pcmci   = PCMCI(dataframe=dataframe, cond_ind_test=parcorr, verbosity=0)
        # results = pcmci.run_pcmci(tau_max=1, pc_alpha=0.05)

        # # evaluate
        # metrics = evaluate_pcmci(true_edges, results, var_names, alpha=0.05)
        # print("PCMCI evaluation:")
        # print(f"  True positives : {metrics['tp']}")
        # print(f"  False positives: {metrics['fp']}")
        # print(f"  False negatives: {metrics['fn']}")
        # print(f"  Precision={metrics['precision']:.3f}, Recall={metrics['recall']:.3f}, F1={metrics['f1']:.3f}")

        # # build predicted edges
        # p_mat = results['p_matrix']
        # N = len(var_names)
        # _, _, tau_max1 = p_mat.shape
        # pred_edges = [
        #     (var_names[i], var_names[j], tau)
        #     for i in range(N) for j in range(N)
        #     for tau in range(1, tau_max1)
        #     if p_mat[j, i, tau] < 0.05
        # ]

        # # plot side by side
        # fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor='white')
        # plot_causal_graph(true_edges,        var_names, "Ground-Truth Graph",       "lightblue",  axes[0])
        # plot_causal_graph(pred_edges,        var_names, "PCMCI-Estimated Graph",     "lightgreen", axes[1])
        # plt.tight_layout()
        # plt.show()
