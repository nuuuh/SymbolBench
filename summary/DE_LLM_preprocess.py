import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipdb
import os
import json
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--eval_type', type=str, default='base', help='Evaluation type')
parser.add_argument('--dataset', type=str, default='strogatz', help='Dataset name')
parser.add_argument('--folder', type=str, default='')
parser.add_argument('--model', type=str, default='Qwen2.5-14B')
parser.add_argument('--modality', type=str, default='textual')
args = parser.parse_args()

eval_type= args.eval_type
dataset = args.dataset

# import ipdb; ipdb.set_trace()

if dataset in 'strogatz_extended':
    gt_file_path = "data/strogatz_extended.json"
    with open("data/strogatz_extended.json", "r") as f:
        gt_data = json.load(f)
    idx_list = list(range(1, 64))
elif dataset == 'physiome':
    with open("data/Physiome/idx_DE", "r") as idx_f:
        idx_list = idx_f.read().strip().split()
    with open("data/Physiome/solved_small_odes.json", "r") as f:
        gt_data = json.load(f)
    idx_list = [int(i) for i in idx_list]

    

if len(args.folder):
    folder = args.folder
else:
    folder = f'runs/{args.model}/{args.modality}/{eval_type}_{dataset}/'

complexity = {}
r2_scores = {}
success = 0

total_files = len(idx_list)
# total_files = 63

results = {'id': [], 'predicted_trees': [], 'duration_fit': [], 'r2': [], 'complexity': []}

# import ipdb; ipdb.set_trace()

for file in gt_data:
    idx = file['id']
    if idx not in idx_list:
        continue
    try:
        data = np.load(folder+f'Expr_{idx}'+'/'+'final_results.npy', allow_pickle=True).item()
        # import ipdb; ipdb.set_trace()
        candidates = data['candidates'][data['candidates']['R2'] >=0.99]
        if candidates.empty:
            candidates = data['candidates'][data['candidates']['R2'] >=0.90]
            if candidates.empty:
                best = data['candidates'].iloc[0]
            else:
                best = candidates.sort_values(by='Complexity', ascending=True).iloc[0]
        else:
            best = candidates.sort_values(by='Complexity', ascending=True).iloc[0]
        # import ipdb; ipdb.set_trace()
        complexity[idx] = best['Complexity']
        r2_scores[idx] = best['R2']
        success += 1

        results['id'].append(idx)
        # import ipdb; ipdb.set_trace()
        # Convert fit_exprs string to list of expression strings
        exprs = [e.strip() for e in str(best['fit_exprs']).strip('[]').split(',')]
        results['predicted_trees'].append('|'.join(exprs))
        results['duration_fit'].append(0)
        results['r2'].append(best['R2'])
        results['complexity'].append(best['Complexity'])

    except:
        print(f"Failed to process {idx}")
        results['id'].append(idx)
        results['predicted_trees'].append(None)
        results['duration_fit'].append(None)
        results['r2'].append(None)
        results['complexity'].append(None)
        continue

sr2 = []
s_complexity = []
for idx, r2 in r2_scores.items():
    if r2>=0:
        sr2.append(r2)
        s_complexity.append(complexity[idx])

sr2 = np.sum(sr2)/total_files
s_complexity = np.sum(s_complexity)/total_files
print(f"Average R2 Score: {sr2}")
print(f"Average Complexity: {s_complexity}")
print(f"Success Rate: {success}/{total_files} = {success/total_files:.2f}")

summary = pd.Series({
    'sr2': sr2,
    's_complexity': s_complexity,
    'success': success,
    'total_files': total_files
    })

# summary.to_csv(f'summary/DE_{eval_type}_summary.csv', index=False)

# import ipdb; ipdb.set_trace()
results_df = pd.DataFrame(results)
results_df.to_csv(f'{folder}/eval_results.csv', index=False)

# import ipdb; ipdb.set_trace()
