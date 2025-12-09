import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

model = "Qwen2.5-14B"

data_type = "CDE"

if data_type == "CDE":
    eval_type = ["naive_strogatz", "base_strogatz", "context_strogatz", "reasoning_strogatz"]
    # eval_type = ["naive_physiome", "base_physiome", "context_physiome", "reasoning_physiome"]
else:
    eval_type = [f"{data_type}_naive", f"{data_type}_base", f"{data_type}_context", f"{data_type}_reasoning"]



converge_epoch = {}

for dataset in eval_type:
    path = f"runs/{model}/textual/{dataset}/"
    converge_epoch[dataset] = []
    for expr in os.listdir(path):
        if expr.startswith("Expr_"):
            try:
                expr_path = os.path.join(path, expr)
                final_results_file = os.path.join(expr_path, "final_results.npy")
                if os.path.exists(final_results_file):
                    data = np.load(final_results_file, allow_pickle=True).item()
                    converge_epoch[dataset].append(data['last_found_at'])
            except Exception as e:
                converge_epoch[dataset].append(0)

for dataset in converge_epoch:
    converge_epoch[dataset] = np.array(converge_epoch[dataset]).mean()

plt.figure(figsize=(10, 6))
plt.bar(converge_epoch.keys(), converge_epoch.values(), color='skyblue')
# plt.xlabel('Dataset Type', fontsize=18)
plt.ylabel('Average Convergence Epoch', fontsize=18)
# plt.title(f'{data_type}', fontsize=18)
plt.xticks(rotation=45, fontsize=16)
plt.yticks(fontsize=16)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(f"summary/{data_type}_convergence.pdf", dpi=300)
plt.savefig(f"summary/{data_type}_convergence.png", dpi=300)
plt.show()

# import ipdb; ipdb.set_trace()
