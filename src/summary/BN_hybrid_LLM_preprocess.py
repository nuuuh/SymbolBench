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
from sympy import simplify, lambdify, Symbol
from scipy.integrate import solve_ivp
from sklearn.metrics import r2_score
from utils.utils import compute_ned
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--eval_type', type=str, default="BN_context", help='Evaluation type') # hybrid_LLM_BN_w_llm BN_context
parser.add_argument('--path', type=str, default="runs/gpt-4o-mini/textual/") # hybrid_output/gpt-4o-mini/textual/ runs/gpt-4o-mini/textual/
args = parser.parse_args()

path = args.path
eval_type = args.eval_type


def final_evaluate(inferred, gt_info):
        train_transitions = gt_info['train_transitions']
        test_transitions = gt_info['test_transitions']

        all_nodes = list(gt_info['train_transitions'][0][0].keys())
        
        # Compute prediction metrics
        id_metrics = test_one_input(train_transitions, all_nodes, inferred)
        ood_metrics = test_one_input(test_transitions, all_nodes, inferred)
        # Compute average complexity of inferred expressions
        # import ipdb; ipdb.set_trace()
        complexities = []
        for node, expr in inferred.items():
            try:
                complexities.append(expr.count_ops())
            except:
                complexities.append(0)
        id_metrics['complexity'] = float(np.mean(complexities)) if complexities else 0.0
        ood_metrics['complexity'] = id_metrics['complexity']

        print(f"ID Metrics: {id_metrics}")
        print(f"OOD Metrics: {ood_metrics}")

        return {'id': id_metrics, 'ood': ood_metrics}

def test_one_input(transitions, all_nodes, inferred):
    """
    Evaluate inferred boolean rules against transitions and compute metrics.
    """
    # import ipdb; ipdb.set_trace()
    TP = FP = FN = TN = 0
    for node, expr in inferred.items():
        sympy_nodes = [Symbol(item) for item in all_nodes]
        f = lambdify(sympy_nodes, expr)
        for x_t, x_tp1 in transitions:
            pred = bool(f(**x_t))
            actual = bool(x_tp1[node])
            if   pred and actual: TP += 1
            elif pred and not actual: FP += 1
            elif not pred and actual: FN += 1
            else: TN += 1
    precision = TP/(TP+FP) if TP+FP>0 else 0
    recall    = TP/(TP+FN) if TP+FN>0 else 0
    f1        = 2*precision*recall/(precision+recall) if precision+recall>0 else 0
    accuracy  = (TP+TN)/(TP+FP+FN+TN) if TP+FP+FN+TN>0 else 0
    bm        = -1 + recall + (TN/(TN+FP) if TN+FP>0 else 0)
    return {'precision':precision,'recall':recall,'f1':f1,'accuracy':accuracy,'BM':bm}


if __name__ == "__main__":
    gt_data_path = "data/BN_small.json"
    with open(gt_data_path, 'r') as f:
        gt_data = json.load(f)
    gt_data = pd.DataFrame(gt_data)

    sample_folder = os.path.join(path, eval_type)
    for folder in os.listdir(sample_folder):
        idx=int(folder.split('_')[-1])
        # import ipdb; ipdb.set_trace()
        if idx not in gt_data['idx'].values:
            continue
        print(f"Processing folder: {folder}")
        try:
            data = np.load(os.path.join(sample_folder, folder, 'final_results.npy'), allow_pickle=True).item()
        except:
            print(f"Skipping {folder} due to missing final_results.npy")
            continue
        
        best = data['candidates'].iloc[0]['fit_expr']
        gt_info = gt_data[gt_data['idx'] == idx].iloc[0]
        metrics = final_evaluate(best, gt_info)
        results = {"exprs": best, "metrics": metrics}
        np.save(os.path.join(sample_folder, folder, 'final_summary.npy'), results)

    