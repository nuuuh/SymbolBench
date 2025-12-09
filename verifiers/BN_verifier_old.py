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

from models.openai_model import OpenAIModel
from models.judge_model import JudgeModel
from .base import BaseVerifier


class BN_Verifier(BaseVerifier):
    def __init__(self, cfg, dataset, logger, device=None, dtype=None):
        super().__init__(cfg)
        self.dataset = dataset
        self.logger = logger
        self.device = device
        self.dtype=dtype
        self.judge_model = None
        if cfg.experiment.symbolic_expression.judge.enabled:
            self.judge_name = cfg.experiment.symbolic_expression.judge.get("name", "gpt-4.1-nano")
            # import ipdb; ipdb.set_trace()
            self.judge_model = JudgeModel(self.judge_name, device, dtype)
            judge_prompt_path = cfg.experiment.symbolic_expression.judge.get("prompt_path", None)
            if judge_prompt_path is not None:
                with open(judge_prompt_path, 'r') as f:
                    self.judge_prompt = json.load(f)
                # import ipdb; ipdb.set_trace()
                self.judge_prompt = '\n '.join(self.judge_prompt.get("judge_prompt"))

        # import ipdb; ipdb.set_trace()

        self.use_context = cfg.experiment.symbolic_expression.get("use_context", False)
        self.reasoning = cfg.experiment.symbolic_expression.get("reasoning", False)

        self.num_vars = dataset.num_vars
        self.num_exprs = dataset.num_eqs
        self.trans_vars = dataset.trans_vars
        self.transitions = dataset.transitions_str
        self.var_list = dataset.var_list
        self.free_vars = dataset.free_vars
        self.domain = dataset.context
        self.var_description = dataset.variable_mapping
        
        self.gt_exprs = dataset.gt_exprs
    
        if cfg.experiment.symbolic_expression.context.image:
            self.image = dataset.image
        else:
            self.image = None
        
        if cfg.experiment.symbolic_expression.additional_prompt == 'image' and self.image is not None:
            self.additional_prompt = self.image
        else:
            self.additional_prompt=None
        
        self.tolerance = cfg.experiment.symbolic_expression.get("tolerance", 0.90)

        # import ipdb; ipdb.set_trace()

        self.candidate_expressions = pd.DataFrame(columns=["expr", "acc", "f1", 'precision', 'recall', 'complexity'])

        self.OOD = cfg.experiment.symbolic_expression.get("OOD", False)

        self.limit_candidates_show=5

        self.timeout_seconds = 15  # adjust as needed

        self.time_series = self.dataset.data['train_transitions']
        self.nodes = self.dataset.var_list
    


    def prepare_prompt(self, base_prompt: str) -> str:
        prompt = '\n '.join(base_prompt)

        prompt = prompt.replace("{transitions}", str(self.transitions))
        prompt = prompt.replace("{trans_vars}", str(self.trans_vars))
        prompt = prompt.replace("{free_vars}", str(self.free_vars))
        prompt = prompt.replace("{num_vars}", str(self.num_vars))
        prompt = prompt.replace("{num_eqs}", str(self.num_exprs))
        prompt = prompt.replace("{var_list}", str(self.var_list))

        if self.use_context:
            context = f"Here is the domain of the problem: {self.domain}.\n"
            context += f"Here are the variables and their descriptions: {self.var_description}.\n"
            self.context = context
            prompt = prompt.replace("{context}", context)

        if not self.candidate_expressions.empty:
            prompt = prompt.replace("{functions}", str(self.candidate_expressions.iloc[:self.limit_candidates_show]))

        # examples: generate a flexible sample Boolean network of dim=num_eqs
        exprs = []
        # determine max terms for each expression (up to 3 or available vars)
        max_terms = min(len(self.var_list) - len(self.free_vars) - 1, 3)
        for var in self.var_list:
            if var in self.free_vars:
                continue
            # choose number of terms in this expression
            n_terms = random.randint(1, max_terms)
            # pick other variables to use
            candidates = [v for v in self.var_list if v != var]
            vars_used = random.sample(candidates, n_terms)
            # build terms, randomly applying NOT
            terms = [f"NOT {v}" if random.random() < 0.5 else v for v in vars_used]
            # combine terms with random AND/OR operators
            expr_body = terms[0]
            for t in terms[1:]:
                op = random.choice(["AND", "OR"])
                expr_body = f"({expr_body} {op} {t})"
            exprs.append(f"{var} <- ({expr_body})")
        eq_str = " ; ".join(exprs)
        if self.reasoning:
            example_str = f'{{"eq": "{eq_str}", "dim": {self.num_exprs}, "reasoning": ... (explain why the generated equations fits the context and transitions)}}'
        else:
            example_str = f'{{"eq": "{eq_str}", "dim": {self.num_exprs}}}'
        prompt = prompt.replace("{example_str}", example_str)

        self.current_prompt = prompt
        # import ipdb; ipdb.set_trace()

        return prompt




    def parse_model_output(self,str_funcs: str):
        import re
        # Extract all dict-like blocks
        dict_blocks = re.findall(r'\{[^}]*\}', str_funcs)
        all_expressions = []
        reasonings = []

        for block in dict_blocks:
            # Extract the equation string after "eq":
            eq_match = re.search(r'"eq"\s*:\s*"([^"]+)"', block)
            if not eq_match:
                self.logger.info(f"Skipping block as no 'eq' found: {block}")
                continue
            eq_str = eq_match.group(1)
            eqs = [e.strip() for e in eq_str.split(";")]
            cleaned_eqs = []
            for eq in eqs:
                # import ipdb; ipdb.set_trace()
                cleaned = self.clean_function(eq)
                if cleaned != "":
                    cleaned_eqs.append(cleaned)
            joined_cleaned = " ; ".join(cleaned_eqs)
            if joined_cleaned:
                try:
                    # import ipdb; ipdb.set_trace()
                    _, function, _ = self.dataset.parse_regulations_to_sympy(joined_cleaned)
                    # import ipdb; ipdb.set_trace()
                    if len(function) != self.num_exprs:
                        self.logger.info(f"Skipping block as number of equations does not match: {block}")
                        continue
                    all_expressions.append(function)
                except Exception as e:
                    self.logger.warning(f"Could not parse line {joined_cleaned}.")
                    self.logger.warning(str(e))
                    continue
            
            reasoning_match = re.search(r'"reasoning"\s*:\s*"([^"]*)"', block)
            # import ipdb; ipdb.set_trace()
            reasonings.append(reasoning_match.group(1) if reasoning_match else None)
        # import ipdb; ipdb.set_trace()
        return all_expressions, reasonings
    
    def check_duplicate(self, expr):
        return False

    
    def update_candidates(self, expr_list, hybrid=False):
        """
        Update the candidate expressions with the new expression.
        """
        # import ipdb; ipdb.set_trace()
        for i, (expr, reasoning) in tqdm(enumerate(tuple(expr_list))):
            try:
                if not self.check_duplicate(expr):
                    # import ipdb; ipdb.set_trace()
                    fit_score = self.fit_score(expr)
                    if fit_score is None:
                        self.logger.info(f"Expression {i} is not valid.")
                        continue
                    # import ipdb; ipdb.set_trace()
                    tmp  ={"expr": str(expr), "reasoning": reasoning}
                    tmp.update(**fit_score)
                    if self.judge_model is not None:
                        tmp = self.score_quality(tmp)
                    self.candidate_expressions = pd.concat([self.candidate_expressions, pd.DataFrame({k: [v] for k, v in tmp.items()})], ignore_index=True)

                    # import ipdb; ipdb.set_trace()
                else:
                    self.logger.info(f"Expression {expr} already exists in the candidate expressions.")
            except:
                self.logger.warning(f"Could not process expression {i}: {expr}.")
                continue
        # import ipdb; ipdb.set_trace()


        self.candidate_expressions = (
            self.candidate_expressions
            .drop_duplicates(subset=['expr'])
            .sort_values(by="f1", ascending=False)
            .reset_index(drop=True)
            .iloc[:100]
        )
        # import ipdb; ipdb.set_trace()
        return self.candidate_expressions
    


    def fit_score(self, exprs):
        # import ipdb; ipdb.set_trace()
        gt_transitions = self.dataset.data['train_transitions']


        TP = FP = FN = TN = 0
        for src, tgt in gt_transitions:
            for node in exprs.keys():
                eq = exprs[node]
                gt_value = tgt[node]
                pred_value = eq.subs(src)
                if pred_value and gt_value:
                    TP += 1
                elif pred_value and not gt_value:
                    FP += 1
                elif not pred_value and gt_value:
                    FN += 1
                else:
                    TN += 1
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        acc = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
        # import ipdb; ipdb.set_trace()

        # Calculate complexity as the number of operators in the expression
        complexity = 0
        for node, eq in exprs.items():
            try:
                complexity += eq.count_ops()
            except Exception as e:
                self.logger.warning(f"Could not count ops for expression {eq}: {e}")
                complexity = None
        # import ipdb; ipdb.set_trace()
        return {"fit_expr": exprs, 
                "acc": acc, 
                "f1": f1,
                "precision": precision,
                "recall": recall,
                "complexity": complexity
                }




    def check_tolerance(self) -> bool:
        if self.candidate_expressions.empty:
            return False
        
        # import ipdb; ipdb.set_trace()

        if self.candidate_expressions.iloc[0]["f1"] >= self.tolerance:
            return True
        else:
            return False

    def score_quality(self, scored_expr):
        # import ipdb; ipdb.set_trace()
        # Prepare and send prompt to judge model
        base_prompt = copy(self.judge_prompt)
        filled_prompt = base_prompt.replace("{transitions}", self.transitions)
        filled_prompt = filled_prompt.replace("{context}", self.context)
        filled_prompt = filled_prompt.replace("{candidate_funcs}", str(scored_expr))

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




    def clean_function(self, function: str) -> str:
        """
        Cleans and validates a list of boolean equations.
        - Parses raw function string into individual equations.
        - Normalizes operators to AND, OR, NOT, XOR.
        - Ensures each expected variable has an equation (excluding free vars).
        """
        # Normalize logical operators and assignment syntax
        function = function.replace("^", "XOR")  # Handle XOR operator
        function = function.replace("AND", "&").replace("OR", "|").replace("NOT ", "~")
        function = function.replace("<-", "=")
        function = re.sub(r"(\w+)\s*->\s*(\w+)", r"(~\1 | \2)", function)
        # Strip list brackets if present
        s = function.strip()
        if s.startswith('[') and s.endswith(']'):
            s = s[1:-1]

        # Split on commas preceding an x<digit> to isolate equations
        parts = re.split(r",\s*(?=x\d+)", s)
        cleaned_eqs = []
        for part in parts:
            eq = part.strip(" '\"")
            # Normalize logical operators
            eq = re.sub(r"\band\b|&&|&", "AND", eq, flags=re.IGNORECASE)
            eq = re.sub(r"\bor\b|\|\||\|", "OR", eq, flags=re.IGNORECASE)
            eq = re.sub(r"\bnot\b|!|~", "NOT", eq, flags=re.IGNORECASE)
            eq = re.sub(r"\bxor\b|\^", "XOR", eq, flags=re.IGNORECASE)
            # Ensure space after NOT operator
            eq = re.sub(r"(NOT)([a-zA-Z])", r"\1 \2", eq)
            # Remove illegal characters
            eq = re.sub(r"[^A-Za-z0-9_=() \t]", "", eq)
            # Balance parentheses
            eq = utils.balance_brackets(eq)
            cleaned_eqs.append(eq.strip())

        # Return cleaned equation list joined by ' ; '
        return ' ; '.join(cleaned_eqs)