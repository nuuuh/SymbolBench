import ipdb
import torch
import numpy as np
import os
import matplotlib.pyplot as plt

exp = 'DE_reasoning'

path = f'runs/Qwen2.5-14B/textual/{exp}' # case_check
index = 61 # 8 27 33 43

R2 = {}
complexity = {}
context_alignment = {}
scientific_plausibility = {}
data = np.load(f'{path}/Expr_{index}/final_results.npy', allow_pickle=True).item()

for idx, item in data['epoch_scores'].items():
    # import ipdb; ipdb.set_trace()
    R2[idx] = item.R2
    complexity[idx] = item.Complexity
    try:
        context_alignment[idx] = item.context_alignment
        scientific_plausibility[idx] = item.scientific_plausibility
    except:
        context_alignment[idx] = 0.0
        scientific_plausibility[idx] = 0.0

# Create a figure with two subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 12))

# Plot R2 in the first subplot
axs[0].plot(list(R2.keys()), list(R2.values()), label='R2', marker='o', color='blue')
axs[0].set_xlabel('Epoch')
axs[0].set_ylabel('R2 Score')
axs[0].set_title('R2 Score Over Epochs')
axs[0].legend()
axs[0].grid()

# Plot other scores in the second subplot
axs[1].plot(list(complexity.keys()), list(complexity.values()), label='Complexity', marker='o', color='orange')
axs[1].plot(list(context_alignment.keys()), list(context_alignment.values()), label='Context Alignment', marker='s', color='green')
axs[1].plot(list(scientific_plausibility.keys()), list(scientific_plausibility.values()), label='Scientific Plausibility', marker='p', color='red')
axs[1].set_xlabel('Epoch')
axs[1].set_ylabel('Score')
axs[1].set_title('Other Scores Over Epochs')
axs[1].legend()
axs[1].grid()

plt.tight_layout()  # Adjust layout to prevent overlapping titles/labels
plt.savefig(f'figures/{exp}_epoch_{index}.png')
plt.show()


# import ipdb; ipdb.set_trace()