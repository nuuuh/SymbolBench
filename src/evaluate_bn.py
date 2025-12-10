import numpy as np
import pygad
import pandas as pd
from datasets import BN_dataset_0
from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from tqdm import tqdm
import os
import sys
import operator
sys.path.append(os.path.join(os.path.dirname(__file__), 'odeformer'))
from odeformer.utils import timeout, MyTimeoutError
import json


import multiprocessing
from deap import creator, base, tools
#from pre_RG import Regulators
import deap
#from HallOfFame import *
#from Gep_simple import gep_simple
from geppy.algorithms.basic import gep_simple
import geppy as gep
from sklearn.cluster import KMeans
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error
from sklearn import preprocessing

import argparse

parser = argparse.ArgumentParser(description='Evaluate BN baseline')
parser.add_argument('--OOD', type=int, default=0, help='Evaluate on OOD data')
args = parser.parse_args()

POP_SIZE = 25
NUM_GENERATIONS = 50
OOD = bool(args.OOD)

def infer_rules_for(target_node, nodes, transitions):
    
    # 1) Build binary input/output lists
    Input_data = [[x_t[var] for var in nodes] for x_t, _ in transitions]
    Out_data   = [x_tp1[target_node] for _, x_tp1 in transitions]

    # 2) Build GEP primitive set
    pset = gep.PrimitiveSet("Main", input_names=nodes)
    pset.add_function(operator.and_, 2)
    pset.add_function(operator.or_, 2)
    pset.add_function(operator.not_, 1)

    # 3) DEAP creator
    creator.create("FitnessMin", base.Fitness, weights=(1, -1))
    creator.create("Individual", gep.Chromosome, fitness=creator.FitnessMin)

    # 4) Toolbox setup
    h = 5
    toolbox = gep.Toolbox()
    toolbox.register("gene_gen", gep.Gene, pset=pset, head_length=h)
    toolbox.register("individual", creator.Individual,
                     gene_gen=toolbox.gene_gen, n_genes=1, linker=None)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("compile", gep.compile_, pset=pset)

    def evaluate(ind):
        func = toolbox.compile(ind)
        # count exact matches
        n_correct = sum(func(*inp) == out for inp, out in zip(Input_data, Out_data))
        # count distinct real regulators used
        vars_used = {v.name for gene in ind for v in gene.kexpression}
        n_regs    = len(vars_used - {"and_","or_","not_"})
        return (n_correct, n_regs)

    toolbox.register("evaluate", evaluate)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # 5) GEP operators
    toolbox.register("mut_uniform", gep.mutate_uniform,
                     pset=pset, ind_pb=2/(2*h+1))
    toolbox.pbs["mut_uniform"] = 0.5
    toolbox.register("mut_invert", gep.invert,    pb=0.5)
    toolbox.register("mut_is_ts", gep.is_transpose, pb=0.5)
    toolbox.register("mut_ris_ts", gep.ris_transpose, pb=0.5)
    toolbox.register("mut_gene_ts", gep.gene_transpose, pb=0.5)
    toolbox.register("cx_1p", gep.crossover_one_point,   pb=0.5)
    toolbox.register("cx_2p", gep.crossover_two_point,   pb=0.5)
    toolbox.register("cx_gene", gep.crossover_gene,      pb=0.5)

    stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    stats.register("avg", np.mean); stats.register("min", np.min); stats.register("max", np.max)

    # 6) Run GEP
    pop = toolbox.population(n=POP_SIZE)
    hof = tools.HallOfFame(5)
    pop2, log = gep_simple(pop, toolbox,
                           n_generations=NUM_GENERATIONS,
                           n_elites=POP_SIZE//5,
                           stats=stats,
                           hall_of_fame=hof,
                           verbose=False)

    # 7) Extract best individual
    best = sorted(pop2,
                  key=lambda ind: (-ind.fitness.values[0], ind.fitness.values[1]))[0]
    # sympy‐expr
    # rule_str = str(gep.simplify(best))
    # expr = parse_expr(rule_str, transformations=standard_transformations)
    # import ipdb; ipdb.set_trace()
    # get the raw sympy expression
    expr = gep.simplify(best)

    return nodes, expr



# def test_one_input(transitions, nodes, raw_solutions):
#     TP = FP = FN = TN = 0
#     # For each target node, apply its truth-table solution to all transitions
#     for node, sol in raw_solutions.items():
#         for x_t, x_tp1 in transitions:
#             idx = 0
#             for i, var in enumerate(nodes):
#                 if x_t[var]:
#                     idx |= (1 << i)
#             pred = bool(sol[idx])
#             actual = x_tp1[node]
#             if pred and actual:
#                 TP += 1
#             elif pred and not actual:
#                 FP += 1
#             elif not pred and actual:
#                 FN += 1
#             else:
#                 TN += 1
#     precision = TP / (TP + FP) if (TP + FP) > 0 else 0
#     recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
#     f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
#     accuracy  = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
#     bm = -1 + recall + TN/(TN+FP) if (TN + FP) > 0 else 0
#     return {'precision': precision, 'recall': recall, 'f1': f1, 'accuracy': accuracy, "BM": bm}


from sympy import lambdify

def test_one_input(transitions, nodes, inferred):
    """
    transitions : list of (x_t, x_tp1) dicts
    nodes       : list of variable names, in order
    inferred    : dict[node_name] = sympy.Expr
    """
    TP = FP = FN = TN = 0
    # import ipdb; ipdb.set_trace()
    # for each node we have a sympy expression
    for node, expr in inferred.items():
        # turn expr into a Python function f(v1, v2, ..., vk) -> {0,1} or bool
        f = lambdify(nodes, expr)

        for x_t, x_tp1 in transitions:
            # build argument vector in the same order as `nodes`
            cur_args = [int(x_t[var]) for var in nodes]
            pred = bool(f(*cur_args))
            actual = bool(x_tp1[node])

            if   pred and actual: TP += 1
            elif pred and not actual: FP += 1
            elif not pred and actual: FN += 1
            else: TN += 1

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall    = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    accuracy  = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
    bm        = -1 + recall + (TN / (TN + FP) if (TN + FP) > 0 else 0)

    return {
        'precision': precision,
        'recall':    recall,
        'f1':        f1,
        'accuracy':  accuracy,
        'BM':        bm
    }


def process_sample(i, bool_data, transitions_list, bn_dataset=None, OOD=False):
    # import ipdb; ipdb.set_trace()
    idx =bool_data[i]['idx']

    train_transitions = transitions_list[transitions_list.idx==idx]['train_transitions'].item()
    test_transitions = transitions_list[transitions_list.idx==idx]['test_transitions'].item()


    # input_obs_list = []
    # transitions = []
    # for j in range(50):
    #     bool_data[i] = bn_dataset.get_data(i)
    #     data = bool_data[i]['sim'].T.to_dict()
    #     input_obs = bool_data[i]['input_observations']
    #     timepoints = sorted(t for t in data if isinstance(t, int))
    #     trans = [(data[t], data[t+1]) for t in timepoints[:-1]]
    #     transitions += [ t for t in trans if t not in transitions ]
    #     input_obs_list.append(input_obs)

    # train_transitions = transitions[:len(transitions)//2]
    # test_transitions = transitions[len(transitions)//2:]
    # import ipdb; ipdb.set_trace()

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
        regs, expr = infer_rules_for(node, nodes, train_transitions)
        inferred[node] = expr
    complexities = [expr.count_ops() for expr in inferred.values()]
    mean_complexity = np.mean(complexities) if complexities else 0
    # import ipdb; ipdb.set_trace()
    if OOD:
        metrics = test_one_input(test_transitions, nodes, inferred)
    else:
        metrics = test_one_input(train_transitions, nodes, inferred)
    metrics['complexity'] = mean_complexity


    metrics['idx'] = i
    # metrics['input_observations'] = input_obs_list
    # metrics['transitions'] = transitions
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


def main():

    with open(f"data/BN.json", 'r') as f:
        transitions_data = pd.DataFrame(json.load(f))

    bool_data = pd.read_csv("data/boolnet.csv")
    bn_dataset = BN_dataset_0(bool_data)
    bool_data = bn_dataset.get_all_data()

    final_result = {"precision": [], "recall": [], "f1": [], "accuracy": [], "complexity": [], "BM": []}
    all_results = []

    for i in tqdm(range(len(bool_data))):
        try:
            metrics = process_sample(i, bool_data, transitions_list=transitions_data, bn_dataset=bn_dataset, OOD=OOD)
        except MyTimeoutError:
            print(f"Sample {i} timed out and was skipped.")
            all_results.append({
                'idx': i,
                # 'input_observations': None,
                # 'transitions': None,
                'regulations': None,
                'mapping': None,
                'data': None,
                'timepoints': None,
                'obs_vars': None,
                'reg_vars': None,
                'all_vars': None,
                'gt_exprs': None,
                'inferred': None
            })
            continue
        except Exception as e:
            print(f"Sample {i} failed with error: {e}")
            all_results.append({
                'idx': i,
                # 'input_observations': None,
                # 'transitions': None,
                'regulations': None,
                'mapping': None,
                'data': None,
                'timepoints': None,
                'obs_vars': None,
                'reg_vars': None,
                'all_vars': None,
                'gt_exprs': None,
                'inferred': None
            })
            continue
        precision = metrics['precision']
        recall = metrics['recall']
        f1 = metrics['f1']
        accuracy = metrics['accuracy']
        mean_complexity = metrics['complexity']
        mean_bm = metrics['BM']
        print(f"Instance {i}: Precision={precision}, Recall={recall}, F1={f1}, Accuracy={accuracy}, BM: {mean_bm}, Complexity={mean_complexity}")
        final_result['precision'].append(precision)
        final_result['recall'].append(recall)
        final_result['f1'].append(f1)
        final_result['accuracy'].append(accuracy)
        final_result['BM'].append(mean_bm)
        final_result['complexity'].append(mean_complexity)
        all_results.append(metrics)

    print("Average Precision:", np.mean(final_result["precision"]))
    print("Average Recall:   ", np.mean(final_result["recall"]))
    print("Average F1:       ", np.mean(final_result["f1"]))
    print("Average Accuracy: ", np.mean(final_result["accuracy"]))
    print("Average BM:       ", np.mean(final_result['BM']))
    print("Average Complexity:", np.mean(final_result["complexity"]))

    out_file = "baseline_results/evaluate_bn_ood.npy" if OOD else "baseline_results/evaluate_bn_id.npy"
    np.save(out_file, {"evaluations": all_results, "scores": final_result})
    print("Finished.")

if __name__ == "__main__":
    main()

