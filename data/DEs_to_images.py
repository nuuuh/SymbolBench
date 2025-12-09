import numpy as np
import torch
import pandas as pd
import os
import json
import matplotlib.pyplot as plt


data_path = "data/gt_DEs.pt"
save_path = 'data/Strogatz_images/'
data = torch.load(data_path)
# import ipdb; ipdb.set_trace()

# data_path = "data/Physiome/solved_small_odes.json"
# save_path = 'data/Physiome_images/'
# with open(data_path, 'r') as f:
#     data = json.load(f)

data = pd.DataFrame(data)

for i, row in data.iterrows():
    # import ipdb; ipdb.set_trace()
    print(f"Processing {i+1}/{len(data)}: {row['id']}")
    time_series = np.array(row['solutions'][0][0]['y']).T
    plt.figure(figsize=(10, 5))
    plt.plot(time_series)
    plt.legend([f"Variable {j+1}" for j in range(time_series.shape[1])])
    plt.savefig(os.path.join(save_path, f"{row['id']}.png"))
    plt.show()
    plt.clf()  # Clear the figure to free memory

# import ipdb; ipdb.set_trace()