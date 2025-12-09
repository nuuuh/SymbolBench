import matplotlib.pyplot as plt
import numpy as np

# Data
metrics = ['Precision', 'Recall', 'Accuracy', 'F1']
bn_wo = [0.612701555, 0.821763819, 0.718109779, 0.696068193]
bn_w  = [0.62924266,  0.807231481, 0.733087571, 0.700440132]

scm_wo = [0.200668279, 0.119905144, 0.119905144, 0.147381203]
scm_w  = [0.223966427, 0.141116303, 0.141116303, 0.168516254]

# Plot setup
x = np.arange(len(metrics))
width = 0.35

# Colors: muted pastel-friendly
colors = {
    'wo': '#f4d4a4',  # dusty rose
    'w':  '#f7e3af'   # pale teal
}

fig, axs = plt.subplots(1, 2, figsize=(10, 4))

# BN subplot
axs[0].bar(x - width/2, bn_wo, width, label='wo', color=colors['wo'])
axs[0].bar(x + width/2, bn_w, width, label='w', color=colors['w'])
axs[0].set_title('BN Qwen3 0.6B', fontsize=18, fontweight='bold')
axs[0].set_xticks(x)
axs[0].set_xticklabels(metrics, fontsize=14)
axs[0].set_ylabel('Score', fontsize=16)
axs[0].legend(title='Thinking', fontsize=14, title_fontsize=10)
axs[0].set_ylim(0.55, 0.85)  # zoom in to highlight small differences
axs[0].grid(axis='y', linestyle=':', alpha=0.5)

# Adding annotations to highlight improvements
for i, (wo, w) in enumerate(zip(bn_wo, bn_w)):
    axs[0].text(i - width/2, wo + 0.005, f'{wo:.2f}', ha='center', fontsize=13, color='black')
    axs[0].text(i + width/2, w + 0.005, f'{w:.2f}', ha='center', fontsize=13, color='black')
    # axs[0].annotate(f'+{(w - wo):.2f}', xy=(i, max(wo, w) + 0.01), ha='center', fontsize=10, color='green')

# SCM subplot
axs[1].bar(x - width/2, scm_wo, width, label='wo', color=colors['wo'])
axs[1].bar(x + width/2, scm_w, width, label='w', color=colors['w'])
axs[1].set_title('SCM Qwen3 1.7B', fontsize=18, fontweight='bold')
axs[1].set_xticks(x)
axs[1].set_xticklabels(metrics, fontsize=14)
axs[1].set_ylim(0.1, 0.25)  # different y-scale to highlight improvements
axs[1].grid(axis='y', linestyle=':', alpha=0.5)

# Adding annotations to highlight improvements
for i, (wo, w) in enumerate(zip(scm_wo, scm_w)):
    axs[1].text(i - width/2, wo + 0.002, f'{wo:.2f}', ha='center', fontsize=13, color='black')
    axs[1].text(i + width/2, w + 0.002, f'{w:.2f}', ha='center', fontsize=13, color='black')
    # axs[1].annotate(f'+{(w - wo):.2f}', xy=(i, max(wo, w) + 0.005), ha='center', fontsize=10, color='green')

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("summary/compute_SCM_BN.pdf", dpi=300)
plt.savefig("summary/compute_SCM_BN.png", dpi=300)
plt.show()