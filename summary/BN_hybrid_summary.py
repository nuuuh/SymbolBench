import numpy as np
import pandas as pd
import os
import json
from copy import deepcopy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--eval_type', type=str, default="BN_context", help='Evaluation type') # hybrid_LLM_BN_w_llm BN_context
parser.add_argument('--path', type=str, default="runs/gpt-4o-mini/textual/") # hybrid_output/gpt-4o-mini/textual/ runs/gpt-4o-mini/textual/
args = parser.parse_args()

path = args.path

if len(args.eval_type) ==0:
    eval_types = ["hybrid_GP_BN_w_llm_0.01", "hybrid_GP_BN_w_llm_0.1", "hybrid_GP_BN_w_llm_0.5", "hybrid_GP_BN_w_llm_1", "hybrid_GP_BN_w_llm_2", "hybrid_LLM_BN_w_llm_1","hybrid_GP_BN_wo_llm" ]
else:
    eval_types = [args.eval_type]

total_results = []

gt_data_path = "data/BN_small.json"
with open(gt_data_path, 'r') as f:
    gt_data = json.load(f)
gt_data = pd.DataFrame(gt_data)


for eval_type in eval_types:
    results = {'id':{'precision':[], 'recall':[], 'f1':[], 'accuracy':[], 'BM':[], 'complexity':[]},
               'ood':{'precision':[], 'recall':[], 'f1':[], 'accuracy':[], 'BM':[], 'complexity':[]}}

    cur_path = f"{path}/{eval_type}/"
    for file in os.listdir(cur_path):
        idx=int(file.split('_')[-1])
        # import ipdb; ipdb.set_trace()
        if idx not in gt_data['idx'].values:
            continue
        print(f"Processing folder: {file}")
        try:
            data = np.load(f"{cur_path}/{file}/final_summary.npy", allow_pickle=True).item()
        except:
            try:
                data = np.load(f"{cur_path}/{file}/final_results.npy", allow_pickle=True).item()   
            except:
                continue
         
        # import ipdb; ipdb.set_trace()
        results['id']['precision'].append(data['metrics']['id']['precision'])
        results['id']['recall'].append(data['metrics']['id']['recall'])
        results['id']['f1'].append(data['metrics']['id']['f1'])
        results['id']['accuracy'].append(data['metrics']['id']['accuracy'])
        results['id']['BM'].append(data['metrics']['id']['BM'])
        results['id']['complexity'].append(data['metrics']['id']['complexity'])
        results['ood']['precision'].append(data['metrics']['ood']['precision'])
        results['ood']['recall'].append(data['metrics']['ood']['recall'])
        results['ood']['f1'].append(data['metrics']['ood']['f1'])
        results['ood']['accuracy'].append(data['metrics']['ood']['accuracy'])
        results['ood']['BM'].append(data['metrics']['ood']['BM'])
        results['ood']['complexity'].append(data['metrics']['ood']['complexity']) 

        tmp = pd.DataFrame(results)
        tmp.to_csv(f"{cur_path}/{file}/summary.csv", index=False)

    if eval_type != "hybrid_BN_wo_llm":
        total_results.append(deepcopy(results))

    results['id']['precision'] = np.mean(results['id']['precision'])
    results['id']['recall'] = np.mean(results['id']['recall'])
    results['id']['f1'] = np.mean(results['id']['f1'])
    results['id']['accuracy'] = np.mean(results['id']['accuracy'])
    results['id']['BM'] = np.mean(results['id']['BM'])
    results['id']['complexity'] = np.mean(results['id']['complexity'])
    results['ood']['precision'] = np.mean(results['ood']['precision'])
    results['ood']['recall'] = np.mean(results['ood']['recall'])
    results['ood']['f1'] = np.mean(results['ood']['f1'])
    results['ood']['accuracy'] = np.mean(results['ood']['accuracy'])
    results['ood']['BM'] = np.mean(results['ood']['BM'])
    results['ood']['complexity'] = np.mean(results['ood']['complexity'])

    results_df = pd.DataFrame(results)
    results_df.to_csv(f"summary/BNs/hybrid_{eval_type}.csv")
    # import ipdb; ipdb.set_trace()


# best_results = {'id':{'precision':[], 'recall':[], 'f1':[], 'accuracy':[], 'BM':[], 'complexity':[]},
#                'ood':{'precision':[], 'recall':[], 'f1':[], 'accuracy':[], 'BM':[], 'complexity':[]}}
# for idx in range(0,29):

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


# best_results['id']['precision'] = np.mean(best_results['id']['precision'])
# best_results['id']['recall'] = np.mean(best_results['id']['recall'])
# best_results['id']['f1'] = np.mean(best_results['id']['f1'])
# best_results['id']['accuracy'] = np.mean(best_results['id']['accuracy'])
# best_results['id']['BM'] = np.mean(best_results['id']['BM'])
# best_results['id']['complexity'] = np.mean(best_results['id']['complexity'])
# best_results['ood']['precision'] = np.mean(best_results['ood']['precision'])
# best_results['ood']['recall'] = np.mean(best_results['ood']['recall'])
# best_results['ood']['f1'] = np.mean(best_results['ood']['f1'])
# best_results['ood']['accuracy'] = np.mean(best_results['ood']['accuracy'])
# best_results['ood']['BM'] = np.mean(best_results['ood']['BM'])
# best_results['ood']['complexity'] = np.mean(best_results['ood']['complexity'])

# best_results_df = pd.DataFrame(best_results)
# best_results_df.to_csv("summary/BNs/hybrid_best_results.csv")



