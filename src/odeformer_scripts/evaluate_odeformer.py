import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import odeformer
from odeformer.model import SymbolicTransformerRegressor
import ipdb;
import numpy as np
from odeformer.metrics import r2_score
from odeformer.odebench.strogatz_equations import equations
from odeformer.odebench.solve_and_plot import config, process_equations, solve_equations
import pandas as pd
import time
import json
import torch


# dataset = "solved_small_odes"
dataset = "strogatz_extended"

# gt_equations = equations # [equations[26]]
# process_equations(gt_equations)
# solve_equations(gt_equations, config)

gt_equations = torch.load('data/gt_DEs.pt')
# with open('data/Physiome/solved_small_odes.json', 'r') as f:
#     gt_equations = json.load(f)


dstr = SymbolicTransformerRegressor(from_pretrained=True)
model_args = {'beam_size':50, 'beam_temperature':0.1}
dstr.set_model_args(model_args)

save_path = 'experiments/odeformer'
# os.mkdir(save_path, exist_ok=True)

config['t_eval'] = np.linspace(0, 10, 150)


results = []
for idx, gt_func in enumerate(gt_equations):
    try:
        gt_eq_str = gt_func['eq']
        t = config['t_eval']
        train_tra = np.array(gt_func['solutions'][0][0]['y']).T
        # import ipdb; ipdb.set_trace()
        start = time.time()
        candidates = dstr.fit(t, train_tra)
        end = time.time()
        duration = end - start
        pred_func = candidates[0][0]
        dstr.print(n_predictions=1)


        results.append({
            'index': gt_func['id'],
            'groundtruth_equation': gt_eq_str,
            'predicted_equation': pred_func,
            "duration_fit": duration,
        })
    except Exception as e:
        print(f"Error processing equation {idx}: {e}")
        results.append({
            'index': gt_func['id'],
            'groundtruth_equation': gt_func['eq'],
            'predicted_equation': None,
            "duration_fit": None,
        })
        continue
import ipdb; ipdb.set_trace()
# Save results to DataFrame
results_df = pd.DataFrame(results)
results_df.to_csv(save_path+f'/eval_{dataset}.csv', index=False)
