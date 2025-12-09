import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipdb
import os
import json
import ast
from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from sympy import symbols, lambdify
import re
from zss import simple_distance
from anytree import NodeMixin


from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from sympy import symbols

def parse_regulations_to_sympy(regulations_str):
    """
    Parses a string of logical regulations into a dictionary of SymPy expressions.

    Args:
        regulations_str (str): A string containing logical regulations in the format:
            'x1 = x6\nx2 = ( ( ( x10 AND NOT x4 ) OR ( x12 AND NOT x4 ) ) OR ( x9 AND NOT x4 ) )'

    Returns:
        dict: A dictionary where keys are variable names (e.g., 'x1') and values are SymPy expressions.
    """
    # Create a dictionary to store the parsed expressions
    sympy_expressions = {}

    # Replace logical operators with Python-compatible syntax
    regulations_str = regulations_str.replace("AND", "&").replace("OR", "|").replace("NOT", "~")

    # Split the string into individual regulations
    regulations = regulations_str.split("\n")

    # Extract variable names and their corresponding expressions
    for regulation in regulations:
        if "=" in regulation:
            var, expr = regulation.split("=", 1)
            var = var.strip()
            expr = expr.strip()

            # Parse the expression into SymPy form
            sympy_expressions[var] = parse_expr(expr, transformations=standard_transformations)

    return sympy_expressions


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




def print_summary(results_df):
    mean_f1 = results_df[~results_df['f1'].isna()]['f1'].mean()
    mean_precision = results_df[~results_df['precision'].isna()]['precision'].mean()
    mean_recall = results_df[~results_df['recall'].isna()]['recall'].mean()
    mean_accuracy = results_df[~results_df['accuracy'].isna()]['accuracy'].mean()
    mean_bm = results_df[~results_df['bm'].isna()]['bm'].mean()
    ACC_05 = np.sum(results_df['f1'] >= 0.5)/65
    ACC_07 = np.sum(results_df['f1'] >= 0.7)/65
    ACC_08 = np.sum(results_df['f1'] >= 0.8)/65
    # import ipdb; ipdb.set_trace()  # Debug
    complexity = results_df['complexity'].mean()
    ned = results_df['ned'].mean()

    summary = {
        'mean_f1': mean_f1,
        'mean_precision': mean_precision,
        'mean_recall': mean_recall,
        'mean_accuracy': mean_accuracy,
        'mean_bm': mean_bm,
        'ACC_05': ACC_05,
        'ACC_07': ACC_07,
        'ACC_08': ACC_08,
        'complexity': complexity,
        'ned': ned
    }

    print(f"Mean F1: {mean_f1}, Mean Precision: {mean_precision}, Mean Recall: {mean_recall}, Mean Accuracy: {mean_accuracy}, BM:{mean_bm}, ACC_05: {ACC_05}, ACC_07: {ACC_07}, ACC_08: {ACC_08}, NED: {ned} Complexity: {complexity}")

    return summary


def evaluate_symbolic_expressions(transitions, all_node_names_list, expressions_dict_sympy):
    """
    Evaluates symbolic expressions against transitions.

    Args:
        transitions: List of tuples (x_t, x_tp1), where x_t and x_tp1 are
                     dicts mapping node names to boolean states.
        all_node_names_list: A list of all node names in the system (used for context, not directly in this version).
        expressions_dict_sympy: A dictionary mapping target_node (str) to its
                                Sympy expression.
    Returns:
        A dictionary with metrics: precision, recall, f1, accuracy.
    """
    TP = FP = FN = TN = 0

    if not transitions:
        return {'precision': 0, 'recall': 0, 'f1': 0, 'accuracy': 0, 'bm': -1, 'error': 'No transitions provided'}

    lambdified_exprs = {}
    for target_node, sympy_expr in expressions_dict_sympy.items():
        if sympy_expr is None:
            continue
        
        args_for_lambdify_func = list(sympy_expr.free_symbols)
        
        try:
            if not args_for_lambdify_func:
                func = lambdify([], sympy_expr, modules="sympy")
            else:
                func = lambdify(args_for_lambdify_func, sympy_expr, modules="sympy")
            lambdified_exprs[target_node] = (func, [str(s) for s in args_for_lambdify_func])
        except Exception as e:
            print(f"Error lambdifying expression for node {target_node}: {sympy_expr}. Error: {e}")
            lambdified_exprs[target_node] = (None, [])

    for x_t, x_tp1 in transitions:
        for target_node, sympy_expr in expressions_dict_sympy.items():
            if sympy_expr is None or target_node not in lambdified_exprs or lambdified_exprs[target_node][0] is None:
                continue

            lmbd_func, arg_names_str = lambdified_exprs[target_node]
            
            call_args = {}
            try:
                for arg_name in arg_names_str:
                    if arg_name not in x_t:
                        # This variable required by the expression is not in the current state's keys.
                        # This could indicate an issue with data consistency or expression variables.
                        # For robustness, one might assume a default value (e.g., False) or skip.
                        # Here, we'll raise to indicate a potential problem.
                        raise KeyError(f"Variable '{arg_name}' for node '{target_node}' (expr: {sympy_expr}) not in current state x_t keys: {list(x_t.keys())}")
                    call_args[arg_name] = x_t[arg_name]
                
                if not arg_names_str: # Constant expression
                    pred_val = lmbd_func()
                else:
                    pred_val = lmbd_func(**call_args)
                pred = bool(pred_val)

            except Exception as e:
                # print(f"Error evaluating expression for {target_node} with state {x_t}: {e}")
                continue # Skip this prediction if evaluation fails

            if target_node not in x_tp1:
                # print(f"Warning: Target node {target_node} not in next state x_tp1. Skipping.")
                continue
            
            actual = x_tp1[target_node]

            if pred and actual:
                TP += 1
            elif pred and not actual:
                FP += 1
            elif not pred and actual:
                FN += 1
            elif not pred and not actual:
                TN += 1
            
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    accuracy  = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
    # BM (Bookmaker Informedness)
    tpr = recall
    tnr = TN / (TN + FP) if (TN + FP) > 0 else 0
    bm = tpr + tnr - 1

    return {'precision': precision, 'recall': recall, 'f1': f1, 'accuracy': accuracy, 'bm': bm}




import argparse
import os

parse = argparse.ArgumentParser()
parse.add_argument('--model', type=str, default='Qwen2.5-14B', help='Model to evaluate')
parse.add_argument('--eval_type', type=str, default='naive', help='Evaluation type')
args = parse.parse_args()


save_path = 'summary/BNs/'
os.makedirs(save_path, exist_ok=True)

with open("data/BN.json", "r") as f:
    gt_data_list = json.load(f)
gt_data = pd.DataFrame(gt_data_list)

eval_type = args.eval_type
model = args.model
folder = f'runs/{model}/textual/BN_{eval_type}/'




id_results_summary = []
ood_results_summary = []
failed = 0

for file_name in os.listdir(folder):
    current_idx = int(file_name.split('_')[-1])
    try:
        result_path = os.path.join(folder, file_name, 'final_results.npy')
        if not os.path.exists(result_path):
            print(f"final_results.npy not found in {os.path.join(folder, file_name)}. Skipping.")
            id_results_summary.append({
            'index': current_idx,
            'precision': None,
            'recall': None,
            'f1': None,
            'accuracy': None,
            'bm': None,
            'complexity': None,
            'ned': None
                })
            ood_results_summary.append({
                'index': current_idx,
                'precision': None,
                'recall': None,
                'f1': None,
                'accuracy': None,
                'bm': None,
                'complexity': None,
                'ned': None
                })
            failed += 1

            continue

        result = np.load(result_path, allow_pickle=True).item()
            
        expr_str_from_file = result['candidates'].iloc[0].expr
        
        # Extract index, be robust if file_name format varies
        
        # import ipdb; ipdb.set_trace()  # Debug
        gt_series = gt_data[gt_data.idx == current_idx]
        gt = gt_series.iloc[0]
        gt_expr = gt.regulations
        
        id_transitions = gt.train_transitions
        ood_transitions = gt.test_transitions
        initial_obs = gt.initial_obs

        all_node_names = sorted(list(initial_obs[0].keys()))
        
        local_symbol_dict = {name: symbols(name) for name in all_node_names}

        expressions_sympy_dict = {}
        def parse_expr_dict(expr_str):
            pattern = r"'([^']+)'\s*:\s*([^,}]+)"
            matches = re.findall(pattern, expr_str)
            return {k: v.strip() for k, v in matches}
        parsed_expr_dict_str_values = parse_expr_dict(expr_str_from_file)

        for node_key, sympy_str_val in parsed_expr_dict_str_values.items():
            expressions_sympy_dict[node_key] = parse_expr(sympy_str_val, local_dict=local_symbol_dict, transformations=standard_transformations)

        complexity = result['candidates'].iloc[0].complexity

        # import ipdb; ipdb.set_trace()  # Debug
        gt_exprs = parse_regulations_to_sympy(gt_expr)

        proximity = 0
        counted=0
        for key in gt_exprs:
            try:
                pred_expr = expressions_sympy_dict[key]
                gt_expr = gt_exprs[key]
                # Convert expressions to trees
                # import ipdb; ipdb.set_trace()
                pred_tree = expr_to_tree(pred_expr)
                gt_tree = expr_to_tree(gt_expr)

                # Compute proximity using tree edit distance
                proximity += simple_distance(pred_tree, gt_tree)
                counted += 1
            except:
                continue

        # ID evaluation
        id_metrics = evaluate_symbolic_expressions(id_transitions, all_node_names, expressions_sympy_dict)
        # OOD evaluation
        ood_metrics = evaluate_symbolic_expressions(ood_transitions, all_node_names, expressions_sympy_dict)
        id_result = {'index': current_idx, "predicted": expr_str_from_file, 'complexity': complexity, 'ned': proximity/counted}
        id_result.update(id_metrics)
        id_results_summary.append(id_result)

        ood_result = {'index': current_idx, "predicted": expr_str_from_file, 'complexity': complexity, 'ned': proximity/counted}
        ood_result.update(ood_metrics)
        # import ipdb; ipdb.set_trace()  # Debug
        ood_results_summary.append(ood_result)

    except Exception as e:
        print(f"An unexpected error occurred while processing file {file_name}: {e}")
        import traceback
        traceback.print_exc()
        id_results_summary.append({
            'index': current_idx,
            'precision': None,
            'recall': None,
            'f1': None,
            'accuracy': None,
            'bm': None,
            'complexity': None,
            'ned': None

        })
        ood_results_summary.append({
            'index': current_idx,
            'precision': None,
            'recall': None,
            'f1': None,
            'accuracy': None,
            'bm': None,
            'complexity': None,
            'ned': None
        })
        failed += 1

# import ipdb; ipdb.set_trace()  # Debug

id_results_df = pd.DataFrame(id_results_summary).sort_values(by='index').reset_index(drop=True)
ood_results_df = pd.DataFrame(ood_results_summary).sort_values(by='index').reset_index(drop=True)

print(f"\nTotal files processed: {len(id_results_df)}")
print(f"Failed to process {failed} files due to errors.")


id_results_df.to_csv(f'{save_path}/{model}_{eval_type}_id_summary.csv', index=False)
ood_results_df.to_csv(f'{save_path}/{model}_{eval_type}_ood_summary.csv', index=False)


print(f"ID Results Summary:")
id_summary = print_summary(id_results_df)
print(f"OOD Results Summary")
ood_summary = print_summary(ood_results_df)

results = {
    'id_summary': id_summary,
    'ood_summary': ood_summary
}
results = pd.DataFrame(results).T

results.to_csv(f'{save_path}/{model}_{eval_type}_summary.csv', index=True)

# import ipdb; ipdb.set_trace()  # Debug
