import numpy as np
import pandas as pd
import os
import json
import ipdb
from zss import simple_distance
from anytree import NodeMixin
from tqdm import tqdm


ood_data = np.load(f"baseline_results/evaluate_bn_ood.npy", allow_pickle=True).item()
ood_data = ood_data['evaluations']

id_data = np.load(f"baseline_results/evaluate_bn_id.npy", allow_pickle=True).item()
id_data = id_data['evaluations']

save_path = "summary/BNs/"
os.makedirs(save_path, exist_ok=True)


class LabeledNode(NodeMixin):
    """
    A custom node class that includes a 'label' attribute for compatibility with zss.
    """
    def __init__(self, label, parent=None):
        self.label = label
        self.parent = parent

def expr_to_tree(expr):
    """
    Converts a logical expression into a tree structure compatible with zss.
    """
    if isinstance(expr, str):
        # Create a labeled node for string expressions
        root = LabeledNode(expr)
        return root
    elif hasattr(expr, 'args'):  # For SymPy-like objects
        root = LabeledNode(type(expr).__name__)
        for arg in expr.args:
            child = expr_to_tree(arg)
            child.parent = root
        return root
    else:
        return LabeledNode(str(expr))


eval_data = {'id': id_data, 'ood': ood_data}
summary_results = {'id': {}, 'ood': {}}

for eval_type, data in eval_data.items():
    results_summary = []
    for item in tqdm(data):
        try:
            idx = item['idx']
            precision = item['precision']
            recall = item['recall']
            f1 = item['f1']
            accuracy = item['accuracy']
            bm = item['BM'] if 'BM' in item else None  # BM (Bookmaker Informedness)
            proximity = 0
            complexity = 0
            for key in item['gt_exprs']:
                pred_expr = item['regulations'][key]
                gt_expr = item['gt_exprs'][key]
                # Convert expressions to trees
                # import ipdb; ipdb.set_trace()
                complexity += pred_expr.count_ops()
                pred_tree = expr_to_tree(pred_expr)
                gt_tree = expr_to_tree(gt_expr)
                
                # Compute proximity using tree edit distance
                proximity += simple_distance(pred_tree, gt_tree)
            results_summary.append({
                'index': idx,
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'accuracy': accuracy,
                'bm': bm,
                'complexity': complexity,
                'ned': proximity / len(item['gt_exprs'])
            })
        except Exception as e:
            print(f"index {item['idx']} failed")
            results_summary.append({
                'index': item['idx'],
                'precision': None,
                'recall': None,
                'f1': None,
                'accuracy': None,
                'bm': None,
                'complexity': None,
                'ned': None
            })

    results_df = pd.DataFrame(results_summary)
    # ipdb.set_trace()

    mean_f1 = results_df[~results_df['f1'].isna()].mean()['f1']
    mean_precision = results_df[~results_df['precision'].isna()].mean()['precision']
    mean_recall = results_df[~results_df['recall'].isna()].mean()['recall']
    mean_accuracy = results_df[~results_df['accuracy'].isna()].mean()['accuracy']
    mean_ned = results_df[~results_df['ned'].isna()].mean()['ned']
    mean_bm = results_df[~results_df['bm'].isna()].mean()['bm'] if 'bm' in results_df.columns else None
    mean_complexity = results_df['complexity'].mean() if 'complexity' in results_df.columns else None

    ACC_05 = np.sum(results_df['f1'] >= 0.5)/65
    ACC_07 = np.sum(results_df['f1'] >= 0.7)/65
    ACC_08 = np.sum(results_df['f1'] >= 0.8)/65

    print(f"Evaluation Type: {eval_type}")
    print(f"Mean F1: {mean_f1}, Mean Precision: {mean_precision}, Mean Recall: {mean_recall}, Mean Accuracy: {mean_accuracy}, Mean BM: {mean_bm},  ACC_0.5: {ACC_05}, ACC_0.7: {ACC_07}, ACC_0.8: {ACC_08}, Mean NED: {mean_ned}")

    summary_results[eval_type] = {
        'mean_f1': mean_f1,
        'mean_precision': mean_precision,
        'mean_recall': mean_recall,
        'mean_accuracy': mean_accuracy,
        'mean_bm': mean_bm,
        'ACC_05': ACC_05,
        'ACC_07': ACC_07,
        'ACC_08': ACC_08,
        'mean_ned': mean_ned,
        "complexity": mean_complexity
    }

    results_df.to_csv(f"{save_path}/BN_baseline_{eval_type}_summary.csv", index=False)

# import ipdb; ipdb.set_trace()  # Debug  
summary_df = pd.DataFrame(summary_results).T
summary_df.to_csv(f"{save_path}/BN_baseline_summary.csv")