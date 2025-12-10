import os
import sys
# os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from tigramite import data_processing as pp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr

import ipdb


def generate_synthetic_var(
    N=5, T=500, L_max=2, graph_type="erdos_renyi",
    edge_prob=0.3, coef_range=(0.5, 1.0), noise_std=0.1, seed=1
):
    rng = np.random.RandomState(seed)

    # 1) Build skeleton directed graph (including possible self-loops)
    if graph_type == "erdos_renyi":
        G = nx.DiGraph()
        for i in range(N):
            for j in range(N):  # allow i == j for self-loops
                if rng.rand() < edge_prob:
                    G.add_edge(i, j)
    elif graph_type == "small_world":
        k = min(N-1, max(2, int(edge_prob*(N-1)) + (int(edge_prob*(N-1))%2)))
        U = nx.watts_strogatz_graph(N, k=k, p=edge_prob, seed=seed)
        G = nx.DiGraph()
        for u, v in U.edges():
            if rng.rand() < 0.5: G.add_edge(u, v)
            else:                G.add_edge(v, u)
        # Add random self-loops
        for n in range(N):
            if rng.rand() < edge_prob:
                G.add_edge(n, n)
    elif graph_type == "scale_free":
        m = max(1, int(edge_prob*(N-1)))
        G = nx.barabasi_albert_graph(N, m=m, seed=seed).to_directed()
        # Add random self-loops
        for n in range(N):
            if rng.rand() < edge_prob:
                G.add_edge(n, n)
    else:
        raise ValueError(f"Unsupported graph_type: {graph_type}")

    # 2) Create weight matrices for each lag
    W = {τ: np.zeros((N, N)) for τ in range(1, L_max+1)}
    for τ in W:
        for i, j in G.edges():
            W[τ][j, i] = rng.uniform(*coef_range)

    # 3) Simulate VAR process
    data = np.zeros((T, N))
    data[:L_max] = rng.randn(L_max, N)
    noise = rng.randn(T, N) * noise_std
    for t in range(L_max, T):
        for τ in range(1, L_max+1):
            data[t] += data[t-τ] @ W[τ].T
        data[t] += noise[t]

    # 4) Extract true edges
    var_names = [f"X{i}" for i in range(N)]
    true_edges = [
        (var_names[i], var_names[j], τ)
        for τ in W for i in range(N) for j in range(N)
        if W[τ][j, i] != 0.0
    ]

    return data, true_edges


def wrap_for_tigramite(data):
    N = data.shape[1]
    var_names = [f"X{i}" for i in range(N)]
    return pp.DataFrame(data, var_names=var_names)


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


def plot_causal_graph(edges, var_names, title, node_color, ax):
    """
    Draws a directed-lag graph on the given Matplotlib Axes.
    edges: list of (src, tgt, lag)
    var_names: list of node names
    """
    G = nx.DiGraph()
    for src, tgt, lag in edges:
        G.add_edge(src, tgt, lag=lag)

    pos = nx.spring_layout(G, seed=42)
    nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=600, ax=ax)
    nx.draw_networkx_labels(G, pos, font_weight='bold', ax=ax)
    nx.draw_networkx_edges(G, pos, arrowstyle='-|>', arrowsize=15, ax=ax)
    edge_labels = {(u, v): f"lag {d['lag']}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,
                                 font_color='gray', ax=ax)
    ax.set_title(title, fontsize=12)
    ax.axis('off')



import numpy as np
import pandas as pd
from data.generate_synthetic_scg import generate_synthetic_var

# Parameter grids
N_list = [3, 4, 5]
L_max_list = [1, 2, 3]
graph_types = ["erdos_renyi", "small_world", "scale_free"]
edge_probs = [0.4, 0.6]
coef_ranges = [(-1, 1)]
noise_stds = [0.01, 0.1]
seeds = [1, 2, 3]

T = 150
T_keep = 100

results = []

for N in N_list:
    for L_max in L_max_list:
        for graph_type in graph_types:
            for edge_prob in edge_probs:
                for coef_range in coef_ranges:
                    for noise_std in noise_stds:
                        for seed in seeds:
                            data, true_edges = generate_synthetic_var(
                                N=N, T=T, L_max=L_max, graph_type=graph_type,
                                edge_prob=edge_prob, coef_range=coef_range,
                                noise_std=noise_std, seed=seed
                            )
                            # Only keep the last 100 time steps
                            data_last = data[-T_keep:]
                            results.append({
                                "N": N,
                                "L_max": L_max,
                                "graph_type": graph_type,
                                "edge_prob": edge_prob,
                                "coef_range": coef_range,
                                "noise_std": noise_std,
                                "seed": seed,
                                "data": data_last,
                                "true_edges": true_edges
                            })

# Example: Save metadata (not the arrays) to CSV
meta = pd.DataFrame([
    {k: v for k, v in r.items() if k not in ("data", "true_edges")}
    for r in results
])
np.save("data/synthetic_sgc.npy", results)









# if __name__ == "__main__":
#     # --- 1) Generate and wrap data ---
#     data, true_edges = generate_synthetic_var(
#         N=3, T=150, L_max=2, graph_type="small_world",
#         edge_prob=0.4, coef_range=(-1, 1.0), noise_std=0.01, seed=1
#     )
#     data = data[50:]
#     dataframe = wrap_for_tigramite(data)
#     var_names = dataframe.var_names

#     # --- 2) Run PCMCI ---
#     parcorr = ParCorr(significance='analytic')
#     pcmci = PCMCI(dataframe=dataframe, cond_ind_test=parcorr, verbosity=0)
#     results = pcmci.run_pcmci(tau_max=2, pc_alpha=None)

#     # --- 3) Evaluate performance ---
#     metrics = evaluate_pcmci(true_edges, results, var_names, alpha=0.05)
#     print("PCMCI evaluation:")
#     print(f"  True positives : {metrics['tp']}")
#     print(f"  False positives: {metrics['fp']}")
#     print(f"  False negatives: {metrics['fn']}")
#     print(f"  Precision={metrics['precision']:.3f}, Recall={metrics['recall']:.3f}, F1={metrics['f1']:.3f}")

#     # --- 4) Plot true vs. estimated graphs side by side ---
#     # Extract predicted edges at alpha
#     p_mat = results['p_matrix']
#     N, _, tau_max1 = p_mat.shape
#     pred_edges = [
#         (var_names[i], var_names[j], τ)
#         for i in range(N) for j in range(N)
#         for τ in range(1, tau_max1)
#         if p_mat[j, i, τ] < 0.05
#     ]

#     # ipdb.set_trace()

#     plt.figure(figsize=(12, 5), facecolor='white')
#     plt.plot(data)
#     plt.savefig("synthetic_data.png", dpi=300)
#     plt.show()