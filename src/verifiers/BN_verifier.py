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

        if hybrid:
            self.expand_candidates()


        return self.candidate_expressions
    

    def expand_candidates(self):
        """
        Use genetic programming to expand and refine the current candidate expressions.
        The initial population for GP is seeded from self.candidate_expressions.
        """
        import operator
        from deap import creator, base, tools
        import geppy as gep
        from geppy.algorithms.basic import gep_simple
        from sympy.parsing.sympy_parser import parse_expr, standard_transformations

        if self.candidate_expressions.empty:
            return

        POP_SIZE = 15
        NUM_GENERATIONS = 15
        HALL_OF_FAME_SIZE = POP_SIZE // 5
        
        # 1. Setup GEP environment
        pset = gep.PrimitiveSet("Main", input_names=self.var_list)
        pset.add_function(operator.and_, 2, name='and_')
        pset.add_function(operator.or_, 2, name='or_')
        pset.add_function(operator.not_, 1, name='not_')

        # Using try-except to handle re-creation of types in interactive sessions
        try:
            creator.create("FitnessMax", base.Fitness, weights=(1.0, -1.0)) # (f1_score, complexity)
            creator.create("Individual", gep.Chromosome, fitness=creator.FitnessMax)
        except Exception:
            pass # types already created

        toolbox = gep.Toolbox()
        h = 5  # head length
        toolbox.register("gene_gen", gep.Gene, pset=pset, head_length=h)
        toolbox.register("individual", creator.Individual, gene_gen=toolbox.gene_gen, n_genes=len(self.trans_vars), linker=None)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("compile", gep.compile_, pset=pset)

        def sympy_to_gep_gene(expr, pset, head_length):
            """
            Converts a sympy expression to a GEP gene using post-order traversal.
            This is a simplified converter and might not cover all edge cases.
            """
            # Mapping from sympy functions to pset function names
            sympy_map = {
                sympy.And: 'and_',
                sympy.Or: 'or_',
                sympy.Not: 'not_',
                # sympy.Xor: 'xor'  # Assuming xor is in pset
            }

            k_expression = []

            def build_k_expr_postorder(sub_expr):
                if isinstance(sub_expr, sympy.Symbol):
                    k_expression.append(pset[sub_expr.name])
                elif isinstance(sub_expr, sympy.logic.boolalg.Boolean):
                    func_name = sympy_map.get(sub_expr.func)
                    if func_name:
                        for arg in sub_expr.args:
                            build_k_expr_postorder(arg)
                        k_expression.append(pset[func_name])
                    else:
                        raise TypeError(f"Unsupported sympy function: {sub_expr.func}")
                else:
                    raise TypeError(f"Unsupported expression type: {type(sub_expr)}")

            build_k_expr_postorder(expr)
            
            # The k-expression is in RPN (postfix). GEP needs prefix.
            # The _simplify_kexpression in geppy expects prefix notation.
            # Let's reverse the RPN to get prefix for functions with arity 2.
            # For a more general solution, a proper tree traversal would be needed.
            # The current build_k_expr was already building a prefix list. Let's re-verify.
            
            prefix_expr = []
            def build_prefix(sub_expr):
                if isinstance(sub_expr, sympy.Symbol):
                    prefix_expr.append(pset[sub_expr.name])
                elif isinstance(sub_expr, sympy.logic.boolalg.Boolean):
                    # handle XOR separately
                    if sub_expr.func is sympy.Xor:
                        # A XOR B is equivalent to (A & ~B) | (~A & B)
                        a, b = sub_expr.args
                        new_expr = sympy.Or(sympy.And(a, sympy.Not(b)), sympy.And(sympy.Not(a), b))
                        # recursively build the prefix for the new expression
                        build_prefix(new_expr)
                        return

                    func_name = sympy_map.get(sub_expr.func)
                    if func_name:
                        prefix_expr.append(pset[func_name])
                        for arg in sub_expr.args:
                            build_prefix(arg)
                    else:
                        raise TypeError(f"Unsupported sympy function: {sub_expr.func}")
                else:
                    raise TypeError(f"Unsupported expression type: {type(sub_expr)}")
            
            build_prefix(expr)

            # Now, create a valid gene array
            gene_array = prefix_expr
            
            # Pad the gene to make it valid
            terminals = [p for p in pset.primitives[0] if p.arity == 0]
            if len(gene_array) > head_length:
                # Truncate if too long, this might lead to invalid expressions but prevents errors
                gene_array = gene_array[:head_length]

            while len(gene_array) < head_length:
                gene_array.append(random.choice(terminals))

            tail_length = head_length * (pset.max_arity - 1) + 1
            while len(gene_array) < head_length + tail_length:
                gene_array.append(random.choice(terminals))
            
            gene = gep.Gene(pset=pset, head_length=head_length)
            gene[:] = gene_array[:head_length + tail_length]
            return gene


        def evaluate(individual):
            # Compile individual to a callable function
            try:
                func = toolbox.compile(individual)
            except (ValueError, TypeError):
                return -1, 1000 # Invalid individual
            
            # Create a dict of sympy expressions from the individual
            exprs = {}
            for i, gene in enumerate(individual):
                target_node = self.trans_vars[i]
                raw_expr = gep.simplify(gene)
                try:
                    sym_expr = parse_expr(str(raw_expr), local_dict={'and_': sympy.And, 'or_': sympy.Or, 'not_': sympy.Not}, transformations=standard_transformations)
                    exprs[target_node] = sym_expr
                except Exception:
                     return -1, 1000 # Failed to parse

            fit_scores = self.fit_score(exprs)
            if fit_scores is None:
                return -1, 1000
            
            return fit_scores['f1'], fit_scores['complexity']

        toolbox.register("evaluate", evaluate)
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.register("mut_uniform", gep.mutate_uniform, pset=pset, ind_pb=1.0/ (2 * h + 1))
        toolbox.register("mut_invert", gep.invert, pb=0.1)
        toolbox.register("mut_is_ts", gep.is_transpose, pb=0.1)
        toolbox.register("mut_ris_ts", gep.ris_transpose, pb=0.1)
        toolbox.register("mut_gene_ts", gep.gene_transpose, pb=0.1)
        toolbox.register("cx_1p", gep.crossover_one_point, pb=0.3)
        toolbox.register("cx_2p", gep.crossover_two_point, pb=0.2)
        toolbox.register("cx_gene", gep.crossover_gene, pb=0.1)

        # 2. Create initial population from self.candidate_expressions
        pop = []
        # Seed population with best candidates so far
        # import  ipdb; ipdb.set_trace()
        for _, row in self.candidate_expressions.iterrows():
            if len(pop) >= POP_SIZE // 2: # Seed with top half
                break
            try:
                # The expression is already a dictionary of sympy expressions
                expr_dict = row['fit_expr']
                genes = []
                for var_name in self.trans_vars:
                    var_symbol = sympy.Symbol(var_name)
                    if var_symbol in expr_dict:
                        gene = sympy_to_gep_gene(expr_dict[var_symbol], pset, h)
                        genes.append(gene)
                    else:
                        # If a var has no equation, create a random gene for it.
                        genes.append(toolbox.gene_gen())
                
                if len(genes) == len(self.trans_vars):
                    # Create an individual and then replace its genes
                    individual = toolbox.individual()
                    for i, gene in enumerate(genes):
                        individual[i] = gene
                    pop.append(individual)

            except Exception as e:
                self.logger.warning(f"Could not convert expression to GEP individual for seeding: {e}")
                continue

        # Fill the rest of the population with random individuals
        remaining_pop_size = POP_SIZE - len(pop)
        if remaining_pop_size > 0:
            pop.extend(toolbox.population(n=remaining_pop_size))

        if not pop:
            self.logger.info("Seeding population from candidates failed. Creating a new random population.")
            pop = toolbox.population(n=POP_SIZE)

        hof = tools.HallOfFame(HALL_OF_FAME_SIZE)

        # 3. Run GEP
        self.logger.info("Expanding candidates using Genetic Programming...")
        pop, log = gep_simple(pop, toolbox, n_generations=NUM_GENERATIONS, n_elites=1, hall_of_fame=hof, verbose=False)
        self.logger.info("Finished GP expansion.")
        # 4. Process results and update self.candidate_expressions
        new_expressions_to_add = []
        for ind in hof:
            exprs = {}
            is_valid = True
            for i, gene in enumerate(ind):
                target_node = self.trans_vars[i]
                try:
                    raw_expr = gep.simplify(gene)
                    sym_expr = parse_expr(str(raw_expr), local_dict={'and_': sympy.And, 'or_': sympy.Or, 'not_': sympy.Not}, transformations=standard_transformations)
                    exprs[target_node] = sym_expr
                except Exception as e:
                    self.logger.warning(f"Could not process gene from HOF: {e}")
                    is_valid = False
                    break
            if is_valid:
                fit_score = self.fit_score(exprs)
                fit_score.update({'expr': str(fit_score['fit_expr']), "reasoning": None})
                new_expressions_to_add.append(fit_score)

        self.candidate_expressions = pd.concat([self.candidate_expressions, pd.DataFrame(new_expressions_to_add)], ignore_index=True)

        self.candidate_expressions = (  self.candidate_expressions
                                        .drop_duplicates(subset=['expr'])
                                        .sort_values(by="f1", ascending=False)
                                        .reset_index(drop=True)
                                        .iloc[:100]
                                        )


    def fit_score(self, exprs):
        # import ipdb; ipdb.set_trace()
        gt_transitions = self.dataset.data['train_transitions']


        TP = FP = FN = TN = 0
        for src, tgt in gt_transitions:
            for node in exprs.keys():
                eq = exprs[node]
                gt_value = tgt[node]
                if isinstance(eq, bool):
                    pred_value = eq
                else:
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
                if isinstance(eq, (bool, sympy.logic.boolalg.BooleanTrue, sympy.logic.boolalg.BooleanFalse)):
                    complexity += 0
                elif isinstance(eq, sympy.Symbol):
                    complexity += 0
                else:
                    complexity += eq.count_ops()
            except Exception as e:
                self.logger.warning(f"Could not count ops for expression {eq}: {e}")
                # Assign a high complexity as a penalty if counting fails
                complexity += 1000
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