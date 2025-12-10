import json
import os
import numpy as np
import matplotlib.pyplot as plt


# CDEs
with open('data/strogatz_extended.json', 'r') as f:
    data1 = json.load(f)
data1 = np.array(data1)

with open('data/Physiome/small_odes.json', 'r') as f:
    data2 = json.load(f)

with open('data/Physiome/idx_DE', 'r') as f:
    data2_idx = f.readlines()
data2_idx = [int(x) for x in data2_idx[0].split(" ")]

dims = {}

for item in data1:
    if item['dim'] not in dims:
        dims[item['dim']] = []
    dims[item['dim']].append(item['dim'])

# import ipdb; ipdb.set_trace()

found = len(data2_idx)
print(f"Found {found} items in data2_idx")
for item in data2:
    if item['id'] not in data2_idx:
        continue
    # print(item['id'])
    found -=1
    if item['dim'] not in dims:
        dims[item['dim']] = []
    dims[item['dim']].append(item['dim'])

# print(found)

plt.figure(figsize=(10, 6))
bar_values = [len(dims[d]) for d in dims.keys()]
bars = plt.bar(dims.keys(), bar_values, color='skyblue')

# Add numbers on top of each bar
for bar, value in zip(bars, bar_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, 
             str(value), ha='center', fontsize=18)

plt.xlabel('Dimension', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.title('Distribution of CDE Dimensions ', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("summary/DE_dims.pdf", dpi=300)
plt.savefig("summary/DE_dims.png", dpi=300)
plt.show()

# import ipdb; ipdb.set_trace()


# --- BN
import pandas as pd
import matplotlib.pyplot as plt
import json

with open('data/BN.json', 'r') as f:
    data = json.load(f)

data = pd.DataFrame(data)
dims = data.num_var.to_list()

# Sort the unique dimensions for consistent ordering
unique_dims = sorted(set(dims))
bar_values = [dims.count(d) for d in unique_dims]

plt.figure(figsize=(10, 6))
bars = plt.bar(unique_dims, bar_values, color='lightcoral')

# Add numbers on top of each bar
for bar, value in zip(bars, bar_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, 
             str(value), ha='center', fontsize=18)

plt.xlabel('Number of Variables', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.title('Distribution of BN Dimensions', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("summary/BN_dims.pdf", dpi=300)
plt.savefig("summary/BN_dims.png", dpi=300)
plt.show()

# import ipdb; ipdb.set_trace()


# --- SCM
with open('data/science_sgc_total.npy', 'rb') as f:
    data = np.load(f, allow_pickle=True)

dims = [item['data'].shape[-1] for item in data if item['data'].shape[-1] < 40]  

# import ipdb; ipdb.set_trace()

# Sort the unique dimensions for consistent ordering
unique_dims = sorted(set(dims))
bar_values = [dims.count(d) for d in unique_dims]

plt.figure(figsize=(16, 12))  # Further increased figure size
bars = plt.bar(unique_dims, bar_values, color='lightgreen')

# Add numbers on top of each bar
for bar, value in zip(bars, bar_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, 
             str(value), ha='center', fontsize=26)  # Further increased font size

plt.xlabel('Number of Variables', fontsize=26)  # Further increased font size
plt.ylabel('Count', fontsize=26)  # Further increased font size
plt.title('Distribution of SCM Dimensions', fontsize=26)  # Further increased font size
plt.xticks(fontsize=22)  # Further increased font size
plt.yticks(fontsize=22)  # Further increased font size
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig("summary/SCM_dims.pdf", dpi=300)
plt.savefig("summary/SCM_dims.png", dpi=300)
plt.show()

# import ipdb; ipdb.set_trace()