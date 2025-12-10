import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipdb
import os

models = ['pcmci', 'lpcmci', 'j-pcmci+']
dataset = 'science' # science synthetic
results = {}
save_path = "summary/SCM/"
os.makedirs(save_path, exist_ok=True)

for model in models:
    failed = 0
    if dataset == 'synthetic':
        base_path = f"baseline_results/evaluate_{model}_{dataset}_sgc.npy"
    else:
        base_path = f"baseline_results/evaluate_{model}_{dataset}_sgc_total.npy"
    data = np.load(base_path, allow_pickle=True).item()
    data = pd.DataFrame(data['evaluations'])

    summary = []
    for i in range(len(data)):
        idx = data.index[i]
        try:
            precision = data.metrics[i]['precision']
            recall = data.metrics[i]['recall']
            f1 = data.metrics[i]['f1']
            # false discovery rate (FDR)
            tp = data.metrics[i]['tp']
            fp = data.metrics[i]['fp']
            fn = data.metrics[i]['fn']
            fdr = fp / (tp + fp) if (tp + fp) > 0 else 0


            true_edges = data.true_edges[i]
            pred_edges = data.pred_edges[i]
            # import ipdb; ipdb.set_trace()

            # Compute SHD
            def edge_list_to_adj_matrix(edges, nodes):
                adj = np.zeros((len(nodes), len(nodes)), dtype=int)
                node_idx = {n: i for i, n in enumerate(nodes)}
                for src, tgt, _ in edges:
                    adj[node_idx[src], node_idx[tgt]] = 1
                return adj

            # Collect all nodes
            all_nodes = set()
            for e in true_edges + pred_edges:
                all_nodes.add(e[0])
                all_nodes.add(e[1])
            all_nodes = sorted(list(all_nodes))

            # Convert to adjacency matrices
            true_adj = edge_list_to_adj_matrix(true_edges, all_nodes)
            pred_adj = edge_list_to_adj_matrix(pred_edges, all_nodes)

            # SHD: number of edge insertions, deletions, or flips needed to transform pred into true (for directed graphs)
            # This is equivalent to the Hamming distance between the adjacency matrices (binary, directed)
            # But ensure both matrices are binary (0/1) and same shape
            diff = np.abs(true_adj - pred_adj)
            shd = np.sum(diff)

            summary.append({
                'index': idx,
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'shd': shd,
                'fdr': fdr,
                'complexity': len(pred_edges)
            })

        except Exception as e:
            failed += 1
            summary.append({
            'index': idx,
            'precision': None,
            'recall': None,
            'f1': None,
            'shd': None,
            'fdr': None,
            'complexity': None
            })
            continue
        # 
    assert len(data) ==190
    summary = pd.DataFrame(summary)
    mean_f1 = summary[~summary['f1'].isna()]['f1'].mean()
    mean_precision = summary[~summary['precision'].isna()]['precision'].mean()
    mean_recall = summary[~summary['recall'].isna()]['recall'].mean()
    mean_shd = summary[~summary['shd'].isna()]['shd'].mean()
    mean_fdr = summary[~summary['fdr'].isna()]['fdr'].mean()
    ACC_05 = np.sum(summary['f1'] >= 0.5)/ len(data)
    ACC_07 = np.sum(summary['f1'] >= 0.7)/ len(data)
    ACC_08 = np.sum(summary['f1'] >= 0.8)/ len(data)
    complexity = summary['complexity'].mean()

    statistics = {
        'mean_f1': mean_f1,
        'mean_precision': mean_precision,
        'mean_recall': mean_recall,
        'mean_shd': mean_shd,
        'mean_fdr': mean_fdr,
        'ACC_05': ACC_05,
        'ACC_07': ACC_07,
        'ACC_08': ACC_08,
        'complexity': complexity
    }

    results[model] = (summary, statistics, failed)



for model in results:
    summary, statistics, failed = results[model]
    print(f"Model: {model}, Failed Instances: {failed}")
    print(f"Mean F1: {statistics['mean_f1']}, Mean Precision: {statistics['mean_precision']}, "
          f"Mean Recall: {statistics['mean_recall']}, Mean SHD: {statistics['mean_shd']}, Mean FDR: {statistics['mean_fdr']}, "
          f"ACC_05: {statistics['ACC_05']}, ACC_07: {statistics['ACC_07']}, "
          f"ACC_08: {statistics['ACC_08']}"
          f", Complexity: {statistics['complexity']}"
          )
    # Save the summary to CSV
    summary.to_csv(f'{save_path}/SCM_{model}_all.csv', index=False)
    # import ipdb; ipdb.set_trace()
    statistics_df = pd.DataFrame(statistics, index=[0])
    statistics_df.to_csv(f'{save_path}/SCM_{model}_summary.csv', index=False)


# import ipdb; ipdb.set_trace()