import numpy as np
import pandas as pd
import ipdb
import os
import argparse

# Helper function to calculate precision, recall, F1 from edge sets
def calculate_prf1(true_edges, pred_edges):
    # Using set of sorted node pairs to treat edges as undirected for P/R/F1
    # Ignores weights for this calculation.
    true_set = set(tuple(sorted(e[:2])) for e in true_edges)
    pred_set = set(tuple(sorted(e[:2])) for e in pred_edges)

    tp = len(true_set.intersection(pred_set))
    fp = len(pred_set - true_set)
    fn = len(true_set - pred_set)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    # false discovery rate (FDR) can also be calculated if needed
    fdr = fp / (tp + fp) if (tp + fp) > 0 else 0

    return precision, recall, f1, fdr

# Helper function for SHD
def edge_list_to_adj_matrix(edges, nodes_list):
    adj = np.zeros((len(nodes_list), len(nodes_list)), dtype=int)
    node_idx_map = {n: i for i, n in enumerate(nodes_list)}
    for src, tgt, _ in edges: # Assuming edges are (src, tgt, weight)
        if src in node_idx_map and tgt in node_idx_map:
            adj[node_idx_map[src], node_idx_map[tgt]] = 1
    return adj




parse = argparse.ArgumentParser()
parse.add_argument('--folder', type=str, default='runs')
parse.add_argument('--model', type=str, default='Qwen2.5-14B', help='Model to evaluate')
parse.add_argument('--eval_type', type=str, default='naive', help='Evaluation type')
args = parse.parse_args()

model = args.model
eval_type = args.eval_type
folder = args.folder


gt_data = np.load("data/science_sgc_total.npy", allow_pickle=True)
gt_data = pd.DataFrame(list(gt_data))

# Determine LLM model name and eval type from folder path for reporting
folder = f'{folder}/{model}/textual/SCM_{eval_type}/' 

# Initialize lists to store metrics for each processed file/instance
summary_results_dict = []
num_failed_instances = 0
num_processed_instances = 0

save_path = 'summary/SCM/'
os.makedirs(save_path, exist_ok=True)


for file_name_or_dir in os.listdir(folder):
    instance_path = os.path.join(folder, file_name_or_dir)
    if not os.path.isdir(instance_path): # Process only directories
        continue

    try:
        # Assuming directory name is the index or contains it
        # Adjust idx extraction if directory naming is different e.g. 'instance_0'
        idx = int(file_name_or_dir.split('_')[-1])

        gt_instance = gt_data[gt_data.idx == idx].iloc[0]
        true_edges_instance = gt_instance.true_edges

        preds_file_path = os.path.join(instance_path, 'final_results.npy')
        if not os.path.exists(preds_file_path):
            # print(f"Warning: final_results.npy not found for {file_name_or_dir}")
            num_failed_instances += 1
            continue

        preds_content = np.load(preds_file_path, allow_pickle=True).item()
        if 'candidates' not in preds_content or preds_content['candidates'].empty:
            # print(f"Warning: No candidates found in {preds_file_path}")
            num_failed_instances += 1
            continue
        
        best_candidate_series = preds_content['candidates'].iloc[0]
        expr_str = best_candidate_series.expr

        exprs_dict_instance = eval(expr_str)
        pred_edges_instance = []
        for src_node, targets_list in exprs_dict_instance.items():
            for tgt_node, weight_val in targets_list:
                src_node_fmt = src_node.replace('x', 'x_') if not src_node.startswith('x_') and 'x' in src_node else src_node
                tgt_node_fmt = tgt_node.replace('x', 'x_') if not tgt_node.startswith('x_') and 'x' in tgt_node else tgt_node
                pred_edges_instance.append((src_node_fmt, tgt_node_fmt, weight_val))

        # Calculate P, R, F1 for the current instance
        p, r, f1, fdr = calculate_prf1(true_edges_instance, pred_edges_instance)

        # Calculate SHD for the current instance
        current_nodes = set()
        for e_src, e_tgt, _ in true_edges_instance + pred_edges_instance:
            current_nodes.add(e_src)
            current_nodes.add(e_tgt)
        
        sorted_nodes = sorted(list(current_nodes))
        
        shd_instance = 0
        if sorted_nodes: # Proceed if there are nodes
            true_adj_instance = edge_list_to_adj_matrix(true_edges_instance, sorted_nodes)
            pred_adj_instance = edge_list_to_adj_matrix(pred_edges_instance, sorted_nodes)
            shd_instance = np.sum(np.abs(true_adj_instance - pred_adj_instance))
        
        num_processed_instances += 1
        # import ipdb; ipdb.set_trace()
        summary_results_dict.append({'index': idx,
                                     'precision': p,
                                     'recall': r,
                                     'f1': f1,
                                     'fdr': fdr,
                                     'shd': shd_instance,
                                     'CI-score':best_candidate_series['CI-score'],
                                     'complexity': best_candidate_series.complexity})


    except Exception as e:
        print(f"Failed to process instance {file_name_or_dir}: {e}")
        num_failed_instances += 1
        summary_results_dict.append({'index': idx,
                                     'precision': None,
                                     'recall': None,
                                     'f1': None,
                                     'fdr': None,
                                     'shd': None,
                                     'CI-score':None,
                                     'complexity': None})
        continue

# import ipdb; ipdb.set_trace()
print(num_failed_instances, " ", num_processed_instances)
# import ipdb; ipdb.set_trace()
results_df = pd.DataFrame(summary_results_dict).sort_values(by='index').reset_index(drop=True)

mean_f1 = results_df['f1'].mean()
mean_precision = results_df['precision'].mean()
mean_recall = results_df['recall'].mean()
mean_shd = results_df['shd'].mean()
mean_fdr = results_df['fdr'].mean()
ACC_05 = np.sum(results_df['f1'] >= 0.5)/len(gt_data) 
ACC_07 = np.sum(results_df['f1'] >= 0.7)/len(gt_data) 
ACC_08 = np.sum(results_df['f1'] >= 0.8)/len(gt_data) 

summary = {
    'mean_f1': mean_f1,
    'mean_precision': mean_precision,
    'mean_recall': mean_recall,
    'mean_shd': mean_shd,
    'mean_fdr': mean_fdr,
    'ACC_05': ACC_05,
    'ACC_07': ACC_07,
    'ACC_08': ACC_08,
    'complexity': results_df['complexity'].mean()
}

# import ipdb; ipdb.set_trace()
print(f"Mean F1: {mean_f1}, Mean Precision: {mean_precision}, Mean Recall: {mean_recall}, Mean SHD: {mean_shd}, Mean FDR: {mean_fdr}, ACC_05: {ACC_05}, ACC_07: {ACC_07}, ACC_08: {ACC_08}")

results_df.to_csv(f'{save_path}/{model}_{eval_type}_all.csv')

summary_df = pd.DataFrame(summary, index=[0])
summary_df.to_csv(f'{save_path}/{model}_{eval_type}_summary.csv', index=False)





