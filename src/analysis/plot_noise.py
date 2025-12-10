import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle

noise_level=(0.01, 0.05, 0.1, 0.2,  0.4)
models=["GPT-4o-mini", "GPT-OSS-20B", "GPT-OSS-120B"] #

SCM_results = {}
DE_results_ID = {}
DE_results_OOD = {}


for model in models:
    SCM_results[model] = []
    DE_results_ID[model] = []
    DE_results_OOD[model] = []
    for noise in noise_level:
        file_path = f"summary/SCM/{model.lower()}_noise_{noise}_all.csv"
        data = pd.read_csv(file_path)
        SCM_results[model].append({"numeric": np.mean(data['f1'].values), "symbolic": np.mean(data['shd'].values)})

        file_path = f"baseline_results/DE_eval_ID/{model.lower()}_textual_strogatz_noise_{noise}.csv"
        data = pd.read_csv(file_path)
        DE_results_ID[model].append({"numeric": np.sum(data[data['R2']>0]['R2'].values)/20, "symbolic": np.mean(data['ned'].values)})

        file_path = f"baseline_results/DE_eval_OOD/{model.lower()}_textual_strogatz_noise_{noise}.csv"
        data = pd.read_csv(file_path)
        DE_results_OOD[model].append({"numeric": np.sum(data[data['R2']>0]['R2'].values)/20, "symbolic": np.mean(data['ned'].values)})

# Define colors and markers for each model
colors = ['#2E86AB', '#A23B72', '#F18F01']  # Blue, Purple, Orange
markers = ['o', 's', '^']  # Circle, Square, Triangle
linestyles = ['-', '--', '-.']  # Solid, Dashed, Dash-dot

# Create subplots for DE and SCM with improved styling
try:
    plt.style.use('seaborn-v0_8-whitegrid')
except:
    try:
        plt.style.use('seaborn-whitegrid')
    except:
        plt.style.use('default')
        
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))

# Plot DE ID results (left subplot)
for i, model in enumerate(models):
    numeric_values = [result["numeric"] for result in DE_results_ID[model]]
    ax1.plot(noise_level, numeric_values, 
             marker=markers[i], 
             color=colors[i],
             label=f"{model}", 
             linestyle=linestyles[i],
             linewidth=2.5,
             markersize=8,
             markerfacecolor=colors[i],
             markeredgecolor='white',
             markeredgewidth=1.5)

ax1.set_ylabel("SR2 (ID)", fontsize=14)
# ax1.set_title("DE ID", fontsize=16, fontweight='bold', pad=20)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_facecolor('#fafafa')
ax1.tick_params(axis='both', which='major', labelsize=12)

# Plot DE OOD results (middle subplot)
for i, model in enumerate(models):
    numeric_values = [result["numeric"] for result in DE_results_OOD[model]]
    ax2.plot(noise_level, numeric_values, 
             marker=markers[i], 
             color=colors[i],
             label=f"{model}", 
             linestyle=linestyles[i],
             linewidth=2.5,
             markersize=8,
             markerfacecolor=colors[i],
             markeredgecolor='white',
             markeredgewidth=1.5)

ax2.set_ylabel("SR2 (OOD)", fontsize=14)
# ax2.set_title("DE OOD", fontsize=16, fontweight='bold', pad=20)
ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_facecolor('#fafafa')
ax2.tick_params(axis='both', which='major', labelsize=12)

# Plot SCM results (right subplot)
for i, model in enumerate(models):
    numeric_values = [result["numeric"] for result in SCM_results[model]]
    ax3.plot(noise_level, numeric_values, 
             marker=markers[i], 
             color=colors[i],
             label=f"{model}", 
             linestyle=linestyles[i],
             linewidth=2.5,
             markersize=8,
             markerfacecolor=colors[i],
             markeredgecolor='white',
             markeredgewidth=1.5)

ax3.set_ylabel("F1", fontsize=14)
# ax3.set_title("SCM", fontsize=16, fontweight='bold', pad=20)
ax3.grid(True, alpha=0.3, linestyle='--')
ax3.set_facecolor('#fafafa')
ax3.tick_params(axis='both', which='major', labelsize=12)

# Add a single legend above both subplots with improved styling
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, 
           loc='upper center', 
           bbox_to_anchor=(0.5, 0.97), 
           ncol=3, 
           fontsize=14,
           frameon=True,
           fancybox=True,
           shadow=True,
           borderpad=1)

plt.tight_layout()
plt.subplots_adjust(top=0.82, wspace=0.25)  # Make room for the legend and bring subplots closer
fig.text(0.5, 0, "Noise Level", ha='center', fontsize=14)
plt.savefig("summary/DE_SCM_noise.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig("summary/DE_SCM_noise.pdf", dpi=300, bbox_inches='tight', facecolor='white')
plt.show()
