import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parse = argparse.ArgumentParser()
parse.add_argument('--model', type=str, default='Qwen2.5-14B', help='Model to evaluate')
parse.add_argument('--eval_type', type=str, default='', help='Evaluation type')
parse.add_argument('--modality', type=str, default='textual')
args = parse.parse_args()

model = args.model
eval_type = args.eval_type
modality = args.modality

try:
    ID_strogatz_path = f"baseline_results/DE_eval_ID/{model}_{modality}_strogatz_{eval_type}.csv"
    OOD_strogatz_path = f"baseline_results/DE_eval_OOD/{model}_{modality}_strogatz_{eval_type}.csv"
    ID_physiome_path = f"baseline_results/DE_eval_ID/{model}_{modality}_physiome_{eval_type}.csv"
    OOD_physiome_path = f"baseline_results/DE_eval_OOD/{model}_{modality}_physiome_{eval_type}.csv"
except:
    try:
        ID_strogatz_path = f"baseline_results/DE_eval_ID/{model}_{modality}_strogatz_extended.csv"
        OOD_strogatz_path = f"baseline_results/DE_eval_OOD/{model}_{modality}_strogatz_extended.csv"
        ID_physiome_path = f"baseline_results/DE_eval_ID/{model}_{modality}_physiome.csv"
        OOD_physiome_path = f"baseline_results/DE_eval_OOD/{model}_{modality}_physiome.csv"
    except:
        raise ValueError("Please provide the correct paths for the evaluation results.")

# try:
#     ID_strogatz_path = f"baseline_results/DE_eval_ID/{model}_strogatz_{eval_type}.csv"
#     OOD_strogatz_path = f"baseline_results/DE_eval_OOD/{model}_strogatz_{eval_type}.csv"
#     ID_physiome_path = f"baseline_results/DE_eval_ID/{model}_physiome_{eval_type}.csv"
#     OOD_physiome_path = f"baseline_results/DE_eval_OOD/{model}_physiome_{eval_type}.csv"
# except:
#     try:
#         ID_strogatz_path = f"baseline_results/DE_eval_ID/{model}_strogatz_extended.csv"
#         OOD_strogatz_path = f"baseline_results/DE_eval_OOD/{model}_strogatz_extended.csv"
#         ID_physiome_path = f"baseline_results/DE_eval_ID/{model}_physiome_extended.csv"
#         OOD_physiome_path = f"baseline_results/DE_eval_OOD/{model}_physiome_extended.csv"
#     except:
#         raise ValueError("Please provide the correct paths for the evaluation results.")


# gt_file_path_1 = f"data/strogatz_extended.json"
# gt_file_path_2 = f"data/Physiome/small_odes.json"

# gt_data_1 = pd.read_json(gt_file_path_1)
# gt_data_1 = pd.DataFrame(gt_data_1)

# gt_data_2 = pd.read_json(gt_file_path_2)
# gt_data_2 = pd.DataFrame(gt_data_2)

# len_gts = len(gt_data_1) + len(gt_data_2)

# dims = gt_data_1['dim'].tolist() + gt_data_2['dim'].tolist()

# plt.figure(figsize=(8, 5))
# counts, bins, patches = plt.hist(dims, bins=range(min(dims), max(dims)+2), edgecolor='black', align='left')
# plt.xlabel('dim')
# plt.ylabel('Count')
# plt.title('Distribution of dims')
# plt.xticks(range(min(dims), max(dims)+1))
# plt.tight_layout()

# # Annotate each bar with its count
# for count, bin_left in zip(counts, bins):
#     if count > 0:
#         plt.text(bin_left, count, str(int(count)), ha='center', va='bottom')

# plt.savefig('summary/DE_dim_distribution.png')
# plt.show()



results = {}
tolerance = 0.9

ID_strogatz_data = pd.read_csv(ID_strogatz_path)
ID_physiome_data = pd.read_csv(ID_physiome_path)
OOD_strogatz_data = pd.read_csv(OOD_strogatz_path)
OOD_physiome_data = pd.read_csv(OOD_physiome_path)

ID_data = pd.concat([ID_strogatz_data, ID_physiome_data], ignore_index=True)
OOD_data = pd.concat([OOD_strogatz_data, OOD_physiome_data], ignore_index=True)



def summarize(data):
    sr2_scores = {1: [], 2: [], 3: [], 4: []}
    pred_com = {1: [], 2: [], 3: [], 4: []}
    gt_com = {1: [], 2: [], 3: [], 4: []}
    symbolic_proximity = {1: [], 2: [], 3: [], 4: []}
    ACC = {1: [], 2: [], 3: [], 4: []}
    # import ipdb; ipdb.set_trace()
    for t in range(1,5):
        eq_t = data[data['num_eqs'] == t]
        ACC[t] = np.mean(eq_t['R2']>= tolerance)
        pred_com[t] = np.mean(eq_t['predicted_complexity'])
        gt_com[t]= np.mean(eq_t['groundtruth_complexity'])
        symbolic_proximity[t] = np.mean(eq_t['ned'])
        num_samples = len(eq_t)
        eq_t = eq_t[eq_t['R2']>=0]
        sr2_scores[t] = np.sum(eq_t['R2'])/num_samples

    result = {
        'sr2': sr2_scores,
        'symbolic_proximity': symbolic_proximity,
        f'ACC_{tolerance}': ACC,
        'pred_com': pred_com,
        'gt_com': gt_com
        }
    return result

ID_results = summarize(ID_data)
OOD_results = summarize(OOD_data)

# Organize ID and OOD summaries into DataFrames
df_ID = pd.DataFrame(ID_results)
df_OOD = pd.DataFrame(OOD_results)
# Rename columns to indicate condition
df_ID = df_ID.rename(columns={col: f"{col}_ID" for col in df_ID.columns})
df_OOD = df_OOD.rename(columns={col: f"{col}_OOD" for col in df_OOD.columns})
# Combine on dimension index
combined_df = pd.concat([df_ID, df_OOD], axis=1)
combined_df.index.name = 'dim'
# Save combined summary
# import ipdb; ipdb.set_trace()
print(combined_df)

combined_df.to_csv(f'summary/DE_summary/DE_{model}_{modality}_{eval_type}_summary.csv')
