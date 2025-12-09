import numpy as np
import json
import os
import pandas as pd
import sympy

def parse_regulations_to_sympy(regulations_str):
        exprs = {}
        variables = []
        free_vars = []
        from sympy.parsing.sympy_parser import parse_expr, standard_transformations
        transformations = standard_transformations
        # import ipdb; ipdb.set_trace()
        if ';' not in regulations_str:
            eqs = regulations_str.strip().split('\n')
        else:
            eqs = regulations_str.strip().split(';')
        for line in eqs:
            if '=' in line:
                lhs, rhs = line.split('=', 1)
                lhs = lhs.strip()
                rhs = rhs.strip()
                rhs = rhs.replace('AND', '&').replace('OR', '|').replace('NOT', '~').replace('!', '~')
                expr = parse_expr(rhs, transformations=transformations)
                exprs[lhs] = expr
                variables.append(lhs)
                free_vars.append(sorted(str(v) for v in expr.free_symbols))
            
        return variables, exprs, free_vars



model_name = "Qwen3-0.6B"
eval_type = "BN_size_compute_w_thinking" #  DE_size_compute_w_thinking size_compute_wo_thinking
path = f"analysis_output/{model_name}/textual/{eval_type}"

with open("data/BN_small.json", 'r') as f:
    gt_data = json.load(f)
gt_data = pd.DataFrame(gt_data)

scores = []
for file in os.listdir(path):
    if file.startswith("Expr"):
        idx = file.split('_')[1]
        result = np.load(f"{path}/{file}/final_results.npy", allow_pickle=True).item()
        total_tokens = np.sum(list(result['epoch_tokens'].values()))
        
        # gt_item = gt_data[gt_data.idx == int(idx)].iloc[0]
        # _, gt_regulations,_ = parse_regulations_to_sympy(gt_item['regulations'])
        # pred_regulations = result['candidates'].iloc[0].fit_expr
        # symbols = [Symbol(f'x{i+1}') for i in range(gt_item['num_var'])]

        try:
            acc = result['candidates'].iloc[0].acc
            f1 = result['candidates'].iloc[0].f1
            precision = result['candidates'].iloc[0].precision
            recall = result['candidates'].iloc[0].recall

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