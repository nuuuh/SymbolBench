import numpy as np
import pandas as pd
import os
from copy import deepcopy
import json

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--eval_type', type=str, default="context_strogatz", help='Evaluation type') # hybrid_LLM_DE_w_llm context_strogatz
parser.add_argument('--path', type=str, default="runs/gpt-4o-mini/textual/") # hybrid_output/gpt-4o-mini/textual/ runs/gpt-4o-mini/textual/
args = parser.parse_args()

path = args.path

if len(args.eval_type) ==0:
    eval_types = ["hybrid_GP_DE_w_llm_0.1","hybrid_GP_DE_wo_llm" ]
else:
    eval_types = [args.eval_type]

gt_data_path = "data/strogatz_extended.json"
with open(gt_data_path, 'r') as f:
    gt_data = json.load(f)
gt_data = pd.DataFrame(gt_data)

total_results = []

for eval_type in eval_types:
    results = {'ID_R2':[], 'OOD_R2':[], 'NED': [], 'complexity':[]}
    cur_path = f"{path}/{eval_type}/"
    for file in os.listdir(cur_path):
        try:
            idx=int(file.split('_')[-1])
            # import ipdb; ipdb.set_trace()
            if idx not in gt_data['id'].values:
                continue
            print(f"Processing folder: {file}")
            try:
                data = np.load(f"{cur_path}/{file}/final_summary.npy", allow_pickle=True).item()
            except:
                try:
                    data = np.load(f"{cur_path}/{file}/final_results.npy", allow_pickle=True).item()   
                except:
                    continue
            results['ID_R2'].append(data['metrics']['ID_R2'])
            results['OOD_R2'].append(data['metrics']['OOD_R2'])
            # import ipdb; ipdb.set_trace()
            results['NED'].append(data['metrics']['ned'])
            results['complexity'].append(data['metrics']['predicted_complexity'])
        except:
            results['ID_R2'].append(np.nan)
            results['OOD_R2'].append(np.nan)
            results['NED'].append(np.nan)
            results['complexity'].append(np.nan)

    if eval_type != "hybrid_DE_wo_llm":
        total_results.append(deepcopy(results))
    
    results_df = pd.DataFrame(results)
    summary_results = {'ID_SR2': None, 'OOD_SR2': None, 'NED': None, 'complexity': None}
    summary_results['ID_SR2'] = np.nansum(results_df['ID_R2'][results_df['ID_R2']> 0]) / len(results_df['ID_R2'])
    summary_results['OOD_SR2'] = np.nansum(results_df['OOD_R2'][results_df['OOD_R2']> 0]) / len(results_df['OOD_R2'])
    summary_results['NED'] = np.nanmean(results_df['NED'])
    summary_results['complexity'] = np.nanmean(results_df['complexity'])

    summary_results = pd.Series(summary_results)
    summary_results.to_csv(f"summary/DE_summary/hybrid_{eval_type}.csv")
    # import ipdb; ipdb.set_trace()


# best_results = {'id':{'R2':[], 'NED': [], 'complexity':[]},
#                'ood':{'R2':[], 'NED': [], 'complexity':[]}}

# for idx in range(1,24):

#     eval_idx=0
#     best_f1 = 0

#     for results in total_results:
#         # import ipdb; ipdb.set_trace()
#         eval_f1 = results['id']['f1'][idx]
#         if eval_f1 > best_f1:
#             best_f1 = eval_f1
#             eval_idx = total_results.index(results)
        
#     best_results['id']['precision'].append(total_results[eval_idx]['id']['precision'][idx])
#     best_results['id']['recall'].append(total_results[eval_idx]['id']['recall'][idx])
#     best_results['id']['f1'].append(total_results[eval_idx]['id']['f1'][idx])
#     best_results['id']['accuracy'].append(total_results[eval_idx]['id']['accuracy'][idx])
#     best_results['id']['BM'].append(total_results[eval_idx]['id']['BM'][idx])
#     best_results['id']['complexity'].append(total_results[eval_idx]['id']['complexity'][idx])
#     best_results['ood']['precision'].append(total_results[eval_idx]['ood']['precision'][idx])
#     best_results['ood']['recall'].append(total_results[eval_idx]['ood']['recall'][idx])
#     best_results['ood']['f1'].append(total_results[eval_idx]['ood']['f1'][idx])
#     best_results['ood']['accuracy'].append(total_results[eval_idx]['ood']['accuracy'][idx])
#     best_results['ood']['BM'].append(total_results[eval_idx]['ood']['BM'][idx])
#     best_results['ood']['complexity'].append(total_results[eval_idx]['ood']['complexity'][idx])


# best_results['id']['precision'] = np.nanmean(best_results['id']['precision'])
# best_results['id']['recall'] = np.nanmean(best_results['id']['recall'])
# best_results['id']['f1'] = np.nanmean(best_results['id']['f1'])
# best_results['id']['accuracy'] = np.nanmean(best_results['id']['accuracy'])
# best_results['id']['BM'] = np.nanmean(best_results['id']['BM'])
# best_results['id']['complexity'] = np.nanmean(best_results['id']['complexity'])
# best_results['ood']['precision'] = np.nanmean(best_results['ood']['precision'])
# best_results['ood']['recall'] = np.nanmean(best_results['ood']['recall'])
# best_results['ood']['f1'] = np.nanmean(best_results['ood']['f1'])
# best_results['ood']['accuracy'] = np.nanmean(best_results['ood']['accuracy'])
# best_results['ood']['BM'] = np.nanmean(best_results['ood']['BM'])
# best_results['ood']['complexity'] = np.nanmean(best_results['ood']['complexity'])

# best_results_df = pd.DataFrame(best_results)
# best_results_df.to_csv("summary/BNs/hybrid_best_results.csv")



