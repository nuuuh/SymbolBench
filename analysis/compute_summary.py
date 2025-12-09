import numpy as np
import json
import os
import pandas as pd

models = ["Qwen3-0.6B", "Qwen3-1.7B", "Qwen3-4B", "Qwen3-8B", "Qwen3-14B", "o4-mini", "DeepSeek-R1-Distill-Llama-70B-free", "DeepSeek-V3", "gpt-oss-20b", "gpt-oss-120b"] # o4-mini DeepSeek-R1-Distill-Llama-70B-free DeepSeek-V3 gpt-oss-20b gpt-oss-120b

eval_types = ["DE_size_compute_w_thinking", "DE_size_compute_wo_thinking"]


results = {}

for model_name in models:
    for eval_type in eval_types:
        try:
            path = f"analysis_output/{model_name}_{eval_type}summary.csv"
            scores = pd.read_csv(path)
        except:
            print(f"Missing {model_name} {eval_type}")
            continue
        scores = scores[scores.idx<=27]
        sr2 = np.sum(scores[scores['R2']>0]['R2'])/27
        acc_9 = np.mean(scores['R2']>=0.9)

        results[f"{model_name}_{eval_type}"] = {
            "mean_R2": sr2,
            "ACC_09": acc_9,
        }


print(pd.DataFrame(results).T)
