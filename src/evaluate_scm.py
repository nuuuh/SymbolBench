import os
import sys
import argparse

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm

from tigramite.data_processing import DataFrame
from tigramite.pcmci import PCMCI
from tigramite.lpcmci import LPCMCI
from tigramite.rpcmci import RPCMCI
from tigramite.jpcmciplus import JPCMCIplus
from tigramite.independence_tests.parcorr import ParCorr


def wrap_for_tigramite(data):
    N = data.shape[1]
    var_names = [f"X{i}" for i in range(N)]
    return DataFrame(data, var_names=var_names)


def evaluate_model(true_edges, results, var_names, alpha=0.05):
    """
    Compare true_edges (list of (src, tgt, lag)) to PCMCI results
    and compute precision, recall, and F1.
    """
    p_mat = results['p_matrix']
    N, _, tau_max1 = p_mat.shape

    pred = set()
    for i in range(N):              # target index
        for j in range(N):          # source index
            for τ in range(1, tau_max1):
                if p_mat[j, i, τ] < alpha:
                    # import ipdb; ipdb.set_trace()
                    if "x_" in true_edges[0][0]:
                        var_names[j] = var_names[j].replace("X", "x_")
                        var_names[i] = var_names[i].replace("X", "x_")
                    pred.add((var_names[j], var_names[i], τ))

    true_set = set(true_edges)
    tp = len(pred & true_set)
    fp = len(pred - true_set)
    fn = len(true_set - pred)

    precision = tp / (tp + fp) if tp + fp else 0.0
    recall    = tp / (tp + fn) if tp + fn else 0.0
    f1        = (2 * precision * recall / (precision + recall)
                 if precision + recall else 0.0)

    return {"pred_edges": list(pred), 'precision': precision, 'recall': recall, 'f1': f1,
            'tp': tp, 'fp': fp, 'fn': fn, 'pred_count': len(pred)}


def plot_causal_graph(edges, var_names, title, node_color, ax):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate Tigramite baselines on synthetic SGC data."
    )
    parser.add_argument('--baseline', type=str, default='pcmci',
                        choices=['pcmci', 'lpcmci', 'rpcmci', 'j-pcmci+'],
                        help='Which Tigramite baseline to run')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Significance level for edge selection')
    parser.add_argument('--dataset', type=str, default="science_sgc", choices=["science_sgc_total", "synthetic_sgc"],)
    args = parser.parse_args()

    # Resolve data path relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, "data", f"{args.dataset}.npy")
    data = np.load(data_path, allow_pickle=True)

    avg_P, avg_R, avg_F1 = [], [], []

    result = []

    for idx, item in tqdm(enumerate(data)):
        time_series = item['data']
        true_edges  = item['true_edges']
        # import ipdb; ipdb.set_trace()

        dataframe = wrap_for_tigramite(time_series)
        var_names = dataframe.var_names
        N = len(var_names)

        cond_test = ParCorr(significance='analytic')
        try:
            import signal
            class TimeoutException(Exception): pass
            def handler(signum, frame):
                raise TimeoutException()
            signal.signal(signal.SIGALRM, handler)
            timeout_seconds = 300  # 5 minutes, adjust as needed
            signal.alarm(timeout_seconds)
            try:
                if args.baseline == 'pcmci':
                    model = PCMCI(dataframe=dataframe, cond_ind_test=cond_test, verbosity=0)
                    results = model.run_pcmci(tau_max=1, pc_alpha=None)

                elif args.baseline == 'lpcmci':
                    model = LPCMCI(dataframe=dataframe, cond_ind_test=cond_test, verbosity=0)
                    results = model.run_lpcmci(
                        link_assumptions=None,
                        tau_min=0,             # minimum lag to consider
                        tau_max=1,             # maximum lag
                        pc_alpha=args.alpha,   # significance level (must be float)
                        # you can also tune other parameters here if needed
                    )

                elif args.baseline == 'rpcmci':
                    model = RPCMCI(dataframe=dataframe, cond_ind_test=cond_test, verbosity=0)
                    # adjust num_regimes/max_transitions as needed
                    results = model.run_rpcmci(tau_max=1, num_regimes=4, max_transitions=4, pc_alpha=None)

                elif args.baseline == 'j-pcmci+':
                    # classify all variables as 'system'; adjust if you have context variables
                    node_classification = {i: "system" for i in range(N)}
                    model = JPCMCIplus(
                        dataframe=dataframe,
                        node_classification=node_classification,
                        cond_ind_test=cond_test,
                        verbosity=0
                    )
                    results = model.run_jpcmciplus(tau_max=1, pc_alpha=None)

                else:
                    raise ValueError(f"Unknown baseline: {args.baseline}")
            finally:
                signal.alarm(0)
        except Exception as e:
            print(f"Error running {args.baseline}: {e}")
            eval_save = {
                'idx': idx,
                'data': time_series,
                'true_edges': true_edges,
                'pred_edges': None,
                'metrics': None
            }

            if "eq_str" in item:
                eval_save['eq'] = item['eq_str']
            # import ipdb; ipdb.set_trace()
            result.append(eval_save)
            continue


        try:
            metrics = evaluate_model(true_edges, results, var_names, alpha=args.alpha)
        except:
            print(f"Error evaluating")
            eval_save = {
                'idx': idx,
                'data': time_series,
                'true_edges': true_edges,
                'pred_edges': None,
                'metrics': None
            }

            if "eq_str" in item:
                eval_save['eq'] = item['eq_str']
            # import ipdb; ipdb.set_trace()
            result.append(eval_save)
            continue


        avg_P.append(metrics['precision'])
        avg_R.append(metrics['recall'])
        avg_F1.append(metrics['f1'])

        eval_save = {
            'idx': idx,
            'data': time_series,
            'true_edges': true_edges,
            'pred_edges': metrics['pred_edges'],
            'metrics': metrics
        }

        if "eq_str" in item:
            eval_save['eq'] = item['eq_str']
        # import ipdb; ipdb.set_trace()
        result.append(eval_save)

    avg_p = np.mean(avg_P)
    avg_r = np.mean(avg_R)
    avg_f1 = np.mean(avg_F1)
    print(f"Baseline: {args.baseline}")
    print("Average Precision: ", avg_p)
    print("Average Recall:    ", avg_r)
    print("Average F1:        ", avg_f1)

    np.save(f"baseline_results/evaluate_{args.baseline}_{args.dataset}.npy", {"evaluations": result, "scores": {"precision": avg_p, "recall": avg_r, "f1": avg_f1}})