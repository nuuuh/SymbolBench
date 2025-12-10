import numpy as np
import json
import os
import pandas as pd
import sympy
import ast

def compute_metrics(expr_dict, true_edges):
    # Normalize true_edges to match expr_dict format
    true_edges_normalized = [(src.replace('x_', 'x'), tgt.replace('x_', 'x'), lag) for src, tgt, lag in true_edges]

    # Flatten expr_dict into a list of edges
    predicted_edges = []
    for target, sources in expr_dict.items():
        for src, lag in sources:
            predicted_edges.append((src, target, lag))

    # Compute true positives
    true_positives = set(predicted_edges) & set(true_edges_normalized)

    # Compute metrics
    tp = len(true_positives)
    fp = len(predicted_edges) - tp
    fn = len(true_edges_normalized) - tp

    precision = tp / len(predicted_edges) if len(predicted_edges) > 0 else 0
    recall = tp / len(true_edges_normalized) if len(true_edges_normalized) > 0 else 0
    f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    accuracy = tp / len(true_edges_normalized) if len(true_edges_normalized) > 0 else 0

    return precision, recall, f1, accuracy


model_name = "Qwen3-1.7B"
eval_type = "SCM_size_compute_w_thinking" #  DE_size_compute_w_thinking size_compute_wo_thinking
path = f"analysis_output/{model_name}/textual/{eval_type}"

gt_data = np.load("data/SCM_science_small.npy", allow_pickle=True)

scores = []
for file in os.listdir(path):
    if file.startswith("Expr"):
        try:
            idx = file.split('_')[1]
            result = np.load(f"{path}/{file}/final_results.npy", allow_pickle=True).item()
            total_tokens = np.sum(list(result['epoch_tokens'].values()))

            for i in range(len(gt_data)):
                if gt_data[i]['idx'] == int(idx):
                    gt_item = gt_data[i]
                    break
            true_edges = gt_item['true_edges']
        
            expr_str = result['candidates'].iloc[0]['expr']
            expr_dict = ast.literal_eval(expr_str)
            precision, recall, f1, acc = compute_metrics(expr_dict, true_edges)
    
        except:
            scores.append({
            "idx": idx,
            "total_tokens": 0,
            "acc": 0,
            "f1": 0,
            "precision": 0,
            "recall": 0
            })
            continue
        scores.append({
            "idx": idx,
            "total_tokens": total_tokens,
            "acc": acc,
            "f1": f1,
            "precision": precision,
            "recall": recall
            })
        

scores = pd.DataFrame(scores)
scores = scores.sort_values(by='idx', ascending=False)

scores.to_csv(f"analysis_output/{model_name}_{eval_type}summary.csv", index=False)