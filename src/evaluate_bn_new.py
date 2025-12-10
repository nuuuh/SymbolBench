#!/usr/bin/env python3
import numpy as np
import pygad
import pandas as pd
from datasets import BN_dataset_0
from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from tqdm import tqdm
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), 'odeformer'))
from odeformer.utils import timeout, MyTimeoutError

POP_SIZE = 1000
NUM_GENERATIONS = 50
OOD = True
FITNESS_THRESHOLD = 0.60  # Stop GA early when this fitness is reached


def make_fitness(target_node, nodes, transitions):
    def fitness(ga_instance, solution, solution_idx):
        # all variables (including itself) as regulators
        regs = list(range(len(nodes)))
        # solution encodes full truth-table bits
        tt_bits = [int(g) for g in solution]
        # count true positives, false positives, false negatives
        TP = FP = FN = TN = 0
        for x_t, x_tp1 in transitions:
            idx = 0
            for i, r in enumerate(regs):
                if x_t[nodes[r]]:
                    idx |= (1 << i)
            pred = bool(tt_bits[idx])
            actual = x_tp1[target_node]
            # import ipdb; ipdb.set_trace()
            if pred and actual:
                TP += 1
            elif pred and not actual:
                FP += 1
            elif not pred and actual:
                FN += 1
            elif not pred and not actual:
                TN += 1
            
        # compute precision, recall, and F1
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        accuracy = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
        return f1
    return fitness


def infer_rules_for(target_node, nodes, transitions):
    # full regulators = all variables
    fitness_func = make_fitness(target_node, nodes, transitions)
    eval_func = make_fitness(target_node, nodes, transitions)

    def on_generation(ga_instance):
        sol, fit, _ = ga_instance.best_solution()
        avg_f = np.mean(ga_instance.last_generation_fitness)
        gen = ga_instance.generations_completed
        # print(f"Gen {gen:3d} → Best={fit:.4f}, Avg={avg_f:.4f}")
        # Early stopping if fitness threshold met: return True to stop GA
        if fit >= FITNESS_THRESHOLD:
            # print(f"Early stopping: reached fitness {fit:.4f} at gen {gen}")
            return True

    k = len(nodes)
    num_genes = 2**k
    gene_space = [[0,1] for _ in range(num_genes)]

    ga = pygad.GA(
        num_generations=NUM_GENERATIONS,
        num_parents_mating=POP_SIZE // 2,
        sol_per_pop=POP_SIZE,
        num_genes=num_genes,
        gene_type=int,
        gene_space=gene_space,
        fitness_func=fitness_func,
        mutation_percent_genes=10,
        on_generation=on_generation,
        suppress_warnings=True
    )
    ga.run()

    solution, fitness, _ = ga.best_solution()
    print("final fitness:", fitness)
    if fitness is None or fitness == 0:
        raise RuntimeError(f"No predictive rule found for {target_node} (fitness={fitness})")

    # decode truth-table and use all variables as regulators
    regs = list(range(len(nodes)))
    tt_bits = [int(g) for g in solution]

    terms = []
    for key, bit in enumerate(tt_bits):
        if bit:
            lits = []
            for i, r in enumerate(regs):
                var = nodes[r]
                # include or negate based on bit in key
                lits.append(var if ((key >> i) & 1) else f"~{var}")
            terms.append("(" + " & ".join(lits) + ")")
    rule_str = " | ".join(terms) if terms else "False"
    expr = parse_expr(rule_str, transformations=standard_transformations)
    return [nodes[r] for r in regs], expr, solution


def test_one_input(transitions, nodes, raw_solutions):
    TP = FP = FN = TN = 0
    # For each target node, apply its truth-table solution to all transitions
    for node, sol in raw_solutions.items():
        for x_t, x_tp1 in transitions:
            idx = 0
            for i, var in enumerate(nodes):
                if x_t[var]:
                    idx |= (1 << i)
            pred = bool(sol[idx])
            actual = x_tp1[node]
            if pred and actual:
                TP += 1
            elif pred and not actual:
                FP += 1
            elif not pred and actual:
                FN += 1
            else:
                TN += 1
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    accuracy  = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
    # bookmaker's measure (bm)
    bm = -1 + recall + TN/(TN+FP) if (TN + FP) > 0 else 0
    return {'precision': precision, 'recall': recall, 'f1': f1, 'accuracy': accuracy, 'BM': bm}


def process_sample(i, bool_data, bn_dataset, OOD):
    raw_solutions = {}
    input_obs_list = []
    transitions = []
    for j in range(50):
        bool_data[i] = bn_dataset.get_data(i)
        data = bool_data[i]['sim'].T.to_dict()
        input_obs = bool_data[i]['input_observations']
        timepoints = sorted(t for t in data if isinstance(t, int))
        trans = [(data[t], data[t+1]) for t in timepoints[:-1]]
        transitions += [ t for t in trans if t not in transitions ]
        input_obs_list.append(input_obs)

    train_transitions = transitions[:len(transitions)//2]
    test_transitions = transitions[len(transitions)//2:]

    obs_vars = bool_data[i]['inputs']
    reg_vars = bool_data[i]['vars']
    all_vars = reg_vars + obs_vars
    gt_exprs = bool_data[i]['regulations']
    mapping = bool_data[i]['mapping']
    data = bool_data[i]['sim'].T.to_dict()
    timepoints = sorted(t for t in data if isinstance(t, int))
    nodes = sorted(data[timepoints[0]].keys())

    inferred = {}
    for node in nodes:
        if node in obs_vars:
            continue
        regs, expr, sol = infer_rules_for(node, nodes, train_transitions)
        inferred[node] = expr
        raw_solutions[node] = sol
    complexities = [expr.count_ops() for expr in inferred.values()]
    mean_complexity = np.mean(complexities) if complexities else 0
    if OOD:
        metrics = test_one_input(test_transitions, nodes, raw_solutions)
    else:
        metrics = test_one_input(train_transitions, nodes, raw_solutions)
    metrics['complexity'] = mean_complexity


    metrics['idx'] = i
    metrics['input_observations'] = input_obs_list
    metrics['transitions'] = transitions
    metrics['regulations'] = inferred
    metrics['mapping'] = mapping
    metrics['data'] = data
    metrics['timepoints'] = timepoints
    metrics['obs_vars'] = obs_vars
    metrics['reg_vars'] = reg_vars
    metrics['all_vars'] = all_vars
    metrics['gt_exprs'] = gt_exprs
    metrics['inferred'] = inferred

    return metrics

# Decorate with timeout (5 minutes = 300 seconds)
process_sample = timeout(seconds=300)(process_sample)


# Worker wrapper to catch timeouts and errors - MOVED TO TOP LEVEL
def _worker(i, bool_data_arg, bn_dataset_arg, OOD_arg):
    try:
        return process_sample(i, bool_data_arg, bn_dataset_arg, OOD_arg)
    except MyTimeoutError:
        print(f"Sample {i} timed out and was skipped.")
    except Exception as e:
        print(f"Sample {i} failed with error: {e}")
    # return empty metrics on failure
    return {'idx': i, 'precision': None, 'recall': None, 'f1': None, 'accuracy': None, 'complexity': None,
            'input_observations': None, 'transitions': None, 'regulations': None, 'mapping': None,
            'data': None, 'timepoints': None, 'obs_vars': None, 'reg_vars': None, 'all_vars': None,
            'gt_exprs': None, 'inferred': None}


def main():
    import multiprocessing
    from functools import partial

    bool_data = pd.read_csv("data/boolnet.csv")
    bn_dataset = BN_dataset_0(bool_data)
    bool_data = bn_dataset.get_all_data()

    final_result = {"precision": [], "recall": [], "f1": [], "accuracy": [], "BM":[], "complexity": []}
    all_results = []

    # Prepare worker function with necessary arguments
    worker_func = partial(_worker, bool_data_arg=bool_data, bn_dataset_arg=bn_dataset, OOD_arg=OOD)

    # Run in parallel across CPUs
    with multiprocessing.Pool() as pool:
        for metrics in tqdm(pool.imap(worker_func, range(len(bool_data))), total=len(bool_data)):
            all_results.append(metrics)
            if metrics['precision'] is not None:
                i = metrics['idx']
                print(f"Instance {i}: Precision={metrics['precision']}, Recall={metrics['recall']},"
                      f" F1={metrics['f1']}, Accuracy={metrics['accuracy']}, Complexity={metrics['complexity']}")

                final_result['precision'].append(metrics['precision'])
                final_result['recall'].append(metrics['recall'])
                final_result['f1'].append(metrics['f1'])
                final_result['accuracy'].append(metrics['accuracy'])
                final_result['BM'].append(metrics['BM'])
                final_result['complexity'].append(metrics['complexity'])

    # Print averages
    print("Average Precision:", np.mean(final_result["precision"]))
    print("Average Recall:   ", np.mean(final_result["recall"]))
    print("Average F1:       ", np.mean(final_result["f1"]))
    print("Average Accuracy: ", np.mean(final_result["accuracy"]))
    print("Average BM:       ", np.mean(final_result["BM"]))


    out_file = "baseline_results/evaluate_bn_ood.npy" if OOD else "baseline_results/evaluate_bn_id.npy"
    # out_file = 'tmp.npy'
    np.save(out_file, {"evaluations": all_results, "scores": final_result})
    print("Finished.")

if __name__ == "__main__":
    main()

