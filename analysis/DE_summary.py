import numpy as np
import json
import os
import pandas as pd

models = ["DeepSeek-R1-Distill-Llama-70B-free", "DeepSeek-V3", "gpt-oss-20b", "gpt-oss-120b"] # o4-mini DeepSeek-R1-Distill-Llama-70B-free DeepSeek-V3 gpt-oss-20b gpt-oss-120b



for model_name in models:
    eval_type = "DE_size_compute_w_thinking" #  DE_size_compute_w_thinking size_compute_wo_thinking
    path = f"analysis_output/{model_name}/textual/{eval_type}"
    if os.path.isdir(path) is False:
        eval_type = "DE_size_compute_wo_thinking"
        path = f"analysis_output/{model_name}/textual/{eval_type}"
    scores = []
    for file in os.listdir(path):
        if file.startswith("Expr"):
            idx = file.split('_')[1]
            try:
                result = np.load(f"{path}/{file}/final_results.npy", allow_pickle=True).item()
                R2 = result['candidates'].iloc[0].R2
                try:
                    total_tokens = np.sum(list(result['epoch_tokens'].values()))
                except:
                    total_tokens = 0    
            except:
                scores.append({
                "idx": idx,
                "total_tokens": 0,
                "R2": 0
                })
                continue
            scores.append({
                "idx": idx,
                "total_tokens": total_tokens,
                "R2": R2
            })
            

    scores = pd.DataFrame(scores)
    scores = scores.sort_values(by='idx', ascending=False)

    scores.to_csv(f"analysis_output/{model_name}_{eval_type}summary.csv", index=False)