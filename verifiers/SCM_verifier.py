import pandas as pd
import numpy as np
import sympy
from utils import utils
import torch
import json
from copy import copy
from tqdm import tqdm
import random
import re
import pingouin as pg

# from models.openai_model import OpenAIModel
from models.judge_model import JudgeModel
from .base import BaseVerifier


class SCM_Verifier(BaseVerifier):
    def __init__(self, cfg, dataset, logger, device=None, dtype=None):
        super().__init__(cfg)
        # import ipdb; ipdb.set_trace()
        
        self.dataset = dataset
        self.logger = logger
        self.device = device
        self.dtype=dtype
        self.judge_model = None
        if cfg.experiment.symbolic_expression.judge.enabled:
            self.judge_name = cfg.experiment.symbolic_expression.judge.get("name", "gpt-4.1-nano")
            # import ipdb; ipdb.set_trace()
            self.judge_model =JudgeModel(self.judge_name, device, dtype)
            judge_prompt_path = cfg.experiment.symbolic_expression.judge.get("prompt_path", None)
            if judge_prompt_path is not None:
                with open(judge_prompt_path, 'r') as f:
                    self.judge_prompt = json.load(f)
                # import ipdb; ipdb.set_trace()
                self.judge_prompt = '\n '.join(self.judge_prompt.get("judge_prompt"))

        self.use_context = cfg.experiment.symbolic_expression.get("use_context", False)
        self.reasoning = cfg.experiment.symbolic_expression.get("reasoning", False)
        self.ts_encoder = cfg.experiment.symbolic_expression.context.get("ts_encoder", False)

        self.num_vars = dataset.num_vars
        self.var_list = dataset.var_list
        self.domain = dataset.domain
        self.var_description = dataset.var_description
        self.gt_edges = dataset.true_edges
        self.gt_graph_str = dataset.graph_str
        self.time_series = dataset.time_series
        self.laggings = tuple(range(1, dataset.lagging+1))
        self.dim = self.num_vars
        
        self.gt_exprs = dataset.eqs_str
    
        if cfg.experiment.symbolic_expression.context.image:
            self.image = dataset.image
        else:
            self.image = None

        if cfg.get("noise_level", 0) > 0:
            # import ipdb; ipdb.set_trace()
            self.time_series = self.inject_noise(self.time_series, cfg.noise_level)
        
        if cfg.experiment.symbolic_expression.additional_prompt == 'image' and self.image is not None:
            self.additional_prompt = self.image
        elif cfg.experiment.symbolic_expression.additional_prompt == 'time_series' and self.time_series is not None:
            self.additional_prompt = self.time_series
        else:
            self.additional_prompt=None
        
        self.tolerance = cfg.experiment.symbolic_expression.get("tolerance", 0.95)


        self.candidate_expressions = pd.DataFrame(columns=["expr", "graph", 'CI-score', 'complexity'])

        self.OOD = False
        self.limit_candidates_show=10
        self.timeout_seconds = 15  # adjust as needed
    
    def inject_noise(self, data, noise_level=0.01):
        """
        Add Gaussian noise to time series data.
        
        Args:
            data: Time series data array
            noise_level: Standard deviation of noise as fraction of signal std
        
        Returns:
            Noisy data
        """
        # Scale noise relative to signal standard deviation for each variable
        signal_std = np.std(data, axis=0, keepdims=True)
        noise = noise_level * signal_std * np.random.randn(*data.shape)
        return data + noise
    

    def prepare_prompt(self, base_prompt):
        # import ipdb; ipdb.set_trace()
        prompt = '\n '.join(base_prompt)

        var_points = ""
        var_list = [str(x).replace("_","") for x in self.var_list]
        for i, var in enumerate(var_list):
            time_series = str(np.round(self.time_series[:,i],3).tolist())
            if self.use_context and self.ts_encoder and self.time_series is not None:
                var_points += f'{var}: <ts><ts/>; '
                # import ipdb; ipdb.set_trace()
            elif self.use_context and self.image is not None:
                var_points += f"{var}: variable{i+1} in the image; "
            elif self.use_context:
                var_points += f"{var}: {time_series};"
        self.var_points = str(var_points)
        prompt = prompt.replace("{var_points}", str(var_points))
        prompt = prompt.replace("{num_vars}", str(self.num_vars))
        prompt = prompt.replace("{dim}", str(self.dim))
        prompt = prompt.replace("{var_list}", str(var_list))
        prompt = prompt.replace("{lagging_list}", str(self.laggings))

        if self.use_context:
            context = ""
            if self.domain is not None:
                context += f"The domain of this problem is called '{self.domain}'.\n"
            if self.var_description is not None:
                context += f"Here are the variables and their descriptions: {self.var_description}.\n"
            if len(context) == 0:
                context = "No additional context is provided."
            self.context = context

            prompt = prompt.replace("{context}", context)

        if not self.candidate_expressions.empty:
            prompt = prompt.replace("{graphs}", str(self.candidate_expressions.drop(columns=['expr']).iloc[:self.limit_candidates_show]))

        # import ipdb; ipdb.set_trace()
        scm_parts = []
        for idx in range(self.num_vars):
            candidates = list(range(self.num_vars))
            num_deps = random.randint(1, min(4, len(candidates)))
            deps = random.sample(candidates, num_deps)
            dep_strs = [f"x{dep}_(lag={random.choice(self.laggings)})" for dep in deps]
            scm_parts.append(f"x{idx} <- " + ", ".join(dep_strs))
        scm_example = " | ".join(scm_parts)
        # build example dictionary and serialize to JSON
        example_dict = {"scm": scm_example, "dim": self.dim}
        if self.reasoning:
            example_dict["reasoning"] = "...(explain why the generated Structured Causal Model fits)"
        example_str = json.dumps(example_dict)
        # import ipdb; ipdb.set_trace()
        prompt = prompt.replace("{example_str}", example_str)

        self.current_prompt = prompt

        return prompt


    def parse_model_output(self, str_funcs):
        # Extract all dict-like blocks
        dict_blocks = re.findall(r'\{[^}]*\}', str_funcs)
        all_expressions = []
        reasonings = []

        for block in dict_blocks:
            try:
                scm_match = re.search(r'"scm"\s*:\s*"([^"]+)"', block)
                if not scm_match:
                    self.logger.info(f"Skipping block as no 'scm' found: {block}")
                    continue
                scm_str = scm_match.group(1)
                # Split by "|" but keep the separator in the string for later joining
                scms = [e.strip() for e in scm_str.split("|")]
                
                scm_dict = {}

                for scm in scms:
                    # import ipdb; ipdb.set_trace()
                    key = scm.split("<-")[0].strip()
                    value = scm.split("<-")[1].strip()
                    # match all
                    matches = re.findall(r'(\w+)_\(lag=(\d+)\)', value)
                    scm_dict[key] = [(var, int(lag)) for var, lag in matches]
                
                # check the keys and laggings
                if len(scm_dict) != self.num_vars:
                    self.logger.info(f"Skipping block as number of variables {len(scm_dict)} does not match expected {self.num_vars}: {block}")
                    continue
                if any(lag not in self.laggings for deps in scm_dict.values() for _, lag in deps):
                    self.logger.info(f"Skipping block as lagging values {self.laggings} do not match expected: {block}")
                
                all_expressions.append(scm_dict)
                
                reasoning_match = re.search(r'"reasoning"\s*:\s*"([^"]*)"', block)
                # import ipdb; ipdb.set_trace()
                reasonings.append(reasoning_match.group(1) if reasoning_match else None)
            except:
                pass


        return all_expressions, reasonings

    
    def check_duplicate(self, expr):
        return False


    def update_candidates(self, expr_list, hybrid=False):
        """
        Update the candidate expressions with the new expression.
        """
        # import ipdb; ipdb.set_trace()
        for expr, reasoning in tqdm(tuple(expr_list)):
            try:
                if not self.check_duplicate(expr):
                    # import ipdb; ipcdb.set_trace()
                    fit_score = self.fit_score(expr)
                    if fit_score is None:
                        self.logger.info(f"Expression {expr} is not valid.")
                        continue
                    # import ipdb; ipdb.set_trace()
                    tmp  ={"expr": str(expr), 'graph': self.dataset.draw_graph(expr) , "reasoning": reasoning}
                    tmp.update(**fit_score)
                    if self.judge_model is not None:
                        tmp = self.score_quality(tmp)
                    self.candidate_expressions = pd.concat([self.candidate_expressions, pd.DataFrame({k: [v] for k, v in tmp.items()})], ignore_index=True)

                    # import ipdb; ipdb.set_trace()
                else:
                    self.logger.info(f"Expression {expr} already exists in the candidate expressions.")
            except:
                self.logger.error(f"Failed to update candidates with expression {expr}.")
                continue


        self.candidate_expressions = (
            self.candidate_expressions
            .drop_duplicates(subset=['graph'])
            .sort_values(by="CI-score", ascending=False)
            .reset_index(drop=True)
            .iloc[:100]
        )
        # import ipdb; ipdb.set_trace()
        return self.candidate_expressions
    
    def fit_score(self, exprs):
        # Compute mean absolute partial correlation for each proposed edge (controls for other parents)
        ts = self.time_series
        T, _ = ts.shape
        scores = []
        for target, parents in exprs.items():
            i = int(target[1:])
            if not parents:
                continue
            # align all series based on maximum lag
            max_lag = max(lag for _, lag in parents)
            length = T - max_lag
            # build DataFrame with target and each parent series
            data = {'target': ts[max_lag:, i]}
            for idx, (src, lag) in enumerate(parents):
                j = int(src[1:])
                col = f'p{idx}'
                data[col] = ts[max_lag - lag:T - lag, j]
            df = pd.DataFrame(data)
            # compute partial/pearson correlation for each parent
            for idx in range(len(parents)):
                xcol = f'p{idx}'
                covars = [f'p{k}' for k in range(len(parents)) if k != idx]
                try:
                    if covars:
                        res = pg.partial_corr(data=df, x=xcol, y='target', covar=covars, method='pearson')
                        r = res.at[0, 'r']
                    else:
                        r = df[xcol].corr(df['target'])
                except Exception:
                    continue
                scores.append(abs(r))
        if not scores:
            return None
        ci_score = float(np.mean(scores))
        complexity = sum(len(p) for p in exprs.values())

        return {'CI-score': ci_score, 'complexity': complexity}


    def score_quality(self, scored_expr):
        # Prepare and send prompt to judge model
        base_prompt = copy(self.judge_prompt)
        filled_prompt = base_prompt.replace("{var_points}", self.var_points)
        filled_prompt = filled_prompt.replace("{context}", self.context)
        # import ipdb; ipdb.set_trace()
        filled_prompt = filled_prompt.replace("{candidate_exprs}", str(scored_expr))

        str_scores = self.judge_model.generate(filled_prompt)

        # Initialize fixed result structure with required keys
        if self.reasoning:
            final = {
                "context_alignment": 0,
                "scientific_plausibility": 0,
                "conciseness_and_clarity": None,
                "logical_coherence": None
            }
        else:
            final = {
                "context_alignment": 0,
                "scientific_plausibility": 0,
            }
        # Parse numeric scores into raw dict (try JSON then regex)
        raw = {}
        try:
            jm = re.search(r'(\{.*?\})', str_scores, re.DOTALL)
            if jm:
                for k, v in json.loads(jm.group(1)).items():
                    iv = int(v)
                    if not (1 <= iv <= 5): raise ValueError
                    raw[k] = iv
            else:
                raise ValueError
        except Exception:
            for key, val_str in re.findall(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*:\s*([0-9]+)', str_scores):
                iv = int(val_str)
                if 1 <= iv <= 5 and key not in raw:
                    raw[key] = iv
        # Map only the predefined keys
        for k in final.keys():
            if k in raw:
                final[k] = raw[k]
        # Merge into scored_expr and return
        # import ipdb; ipdb.set_trace()
        scored_expr.update(final)
        return scored_expr

    def check_tolerance(self,):
        """
        Check if the candidate expressions are within the tolerance.
        """
        if self.candidate_expressions.empty:
            return False
        # import ipdb; ipdb.set_trace()

        if self.candidate_expressions.iloc[0]["CI-score"] >= self.tolerance:
            return True
        else:
            return False