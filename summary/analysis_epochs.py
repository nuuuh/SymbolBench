import ipdb
import torch
import numpy as np
import os
import matplotlib.pyplot as plt

fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(14, 10))  # Two subplots vertically

# ----------- CDEs
exp_cde = 'reasoning_strogatz'
path_cde = f'runs/ChatTS-14B/textual/{exp_cde}'
index_cde = 27

R2 = {}
complexity_cde = {}
context_alignment_cde = {}
scientific_plausibility_cde = {}
data_cde = np.load(f'{path_cde}/Expr_{index_cde}/final_results.npy', allow_pickle=True).item()

for idx, item in data_cde['epoch_scores'].items():
    R2[idx] = item.R2
    complexity_cde[idx] = item.Complexity
    try:
        context_alignment_cde[idx] = item.context_alignment
        scientific_plausibility_cde[idx] = item.scientific_plausibility
    except:
        context_alignment_cde[idx] = 0.0
        scientific_plausibility_cde[idx] = 0.0

epochs_cde = list(R2.keys())

# Left Y-axis for CDEs
ax1.plot(epochs_cde, list(complexity_cde.values()), label='Complexity', marker='o', markersize=6, color='orange')
ax1.set_ylabel('Complexity', fontsize=24)
ax1.tick_params(axis='y', labelsize=24)
ax1.tick_params(axis='x', labelsize=24)

# Right Y-axis for R2
ax2 = ax1.twinx()
ax2.plot(epochs_cde, list(R2.values()),label=r'$R^2$', marker='D', markersize=6, color='blue')
ax2.set_ylabel(ylabel=r'$R^2$', fontsize=24)
ax2.tick_params(axis='y', labelsize=24)

# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=24, loc='upper left')

# Add internal title for CDE
ax1.text(0.5, 0.9, 'CDE', transform=ax1.transAxes,
         fontsize=22, fontweight='bold', ha='center', va='center')

# ----------- SCMs
exp_scm = 'SCM_reasoning'
path_scm = f'runs/ChatTS-14B/textual/{exp_scm}'
index_scm = 90

ci = {}
complexity_scm = {}
context_alignment_scm = {}
scientific_plausibility_scm = {}
data_scm = np.load(f'{path_scm}/Expr_{index_scm}/final_results.npy', allow_pickle=True).item()

for idx, item in data_scm['epoch_scores'].items():
    ci[idx] = item['CI-score']
    complexity_scm[idx] = item.complexity
    try:
        context_alignment_scm[idx] = item.context_alignment
        scientific_plausibility_scm[idx] = item.scientific_plausibility
    except:
        context_alignment_scm[idx] = 0.0
        scientific_plausibility_scm[idx] = 0.0

epochs_scm = list(ci.keys())

# Left Y-axis for SCMs
ax3.plot(epochs_scm, list(complexity_scm.values()), label='Complexity', marker='o', markersize=6, color='orange')
ax3.set_xlabel('Epoch', fontsize=24)
ax3.set_ylabel('Complexity', fontsize=24)
ax3.tick_params(axis='y', labelsize=24)
ax3.tick_params(axis='x', labelsize=24)

# Right Y-axis for CI-score
ax4 = ax3.twinx()
ax4.plot(epochs_scm, list(ci.values()), label='CI-score', marker='D', markersize=6, color='blue')
ax4.set_ylabel('CI Score', fontsize=24)
ax4.tick_params(axis='y', labelsize=24)
ax4.set_ylim(0, 1.05)

# Combine legends
lines3, labels3 = ax3.get_legend_handles_labels()
lines4, labels4 = ax4.get_legend_handles_labels()
ax3.legend(lines3 + lines4, labels3 + labels4, fontsize=24, loc='upper left', bbox_to_anchor=(0, 0.7))

# Add internal title for SCM
ax3.text(0.5, 0.9, 'SCM', transform=ax3.transAxes,
         fontsize=24, fontweight='bold', ha='center', va='center')

plt.tight_layout(pad=2.0)
plt.savefig('summary/epoch_scores.pdf', dpi=300)
plt.savefig('summary/epoch_scores.png', dpi=300)
plt.show()