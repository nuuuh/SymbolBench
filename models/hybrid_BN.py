import os
import numpy as np
import pandas as pd
import operator
from deap import creator, base, tools
import geppy as gep
from geppy.algorithms.basic import gep_simple
from tqdm import tqdm
from sympy.parsing.sympy_parser import parse_expr, standard_transformations
from sympy import simplify, lambdify
import json
from utils import utils
from transformers import AutoTokenizer, AutoModelForCausalLM
import re
from sympy import Symbol
import random
import torch

from .hf_model import HuggingFaceModel
from .openai_model import OpenAIModel



class Hybrid_BN(object):
    def __init__(self, model_name, device, dtype, cache_dir=None, **kwargs):
        self.set_seed(15)
        self.model_name = model_name
        self.device = device
        self.dtype = dtype
        token = os.environ.get("HF_TOKEN", None)

        if "gpt" not in model_name:
            self.llm = HuggingFaceModel(model_name, device, dtype, cache_dir, **kwargs)
        elif "gpt" in model_name:
            self.llm = OpenAIModel(model_name, device, dtype, cache_dir, **kwargs)


        # # initialize LLM model and tokenizer as before
        # self.tokenizer = AutoTokenizer.from_pretrained(model_name, torch_dtype=dtype, cache_dir=cache_dir, token=token)
        # self.model = AutoModelForCausalLM.from_pretrained(model_name, torch_dtype=dtype, cache_dir=cache_dir, token=token, device_map='auto')
        # if "tokenizer_pad" in kwargs:
        #     self.tokenizer.pad_token = kwargs["tokenizer_pad"]
        # if "tokenizer_padding_side" in kwargs:
        #     self.tokenizer.padding_side = kwargs["tokenizer_padding_side"]

        # self.model.eval()

        # self.temperature = kwargs.get("temperature", 1.0)
        # self.top_k = kwargs.get("top_k", 50)
        # self.top_p = kwargs.get("top_p", 0.9)
        # self.num_beams = kwargs.get("num_beams", 1)
        # self.num_return_sequences = kwargs.get("num_return_sequences", 1)
        # self.max_new_tokens = kwargs.get("max_new_tokens", 256)
        # self.min_new_tokens = kwargs.get("min_new_tokens", 0)
    
        # Boolean network inference setup
        self.POP_SIZE = 15
        self.NUM_GENERATIONS = 15
        self.hall_of_fame = kwargs.get('hof_size', self.POP_SIZE//5)
        # DEAP & GEP primitive set functions
        self.pset = None
        self.toolbox = None
        self.nodes = []
        self.inferred = {}
        self.use_llm=False
        self.llm_weight=1

    def set_seed(self, seed):
        """Set random seed for reproducibility."""
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    # def generate(self, prompt, additional_prompt=None, return_prompt=False, temperature=None, max_new_tokens=None):
    #     if temperature is None:
    #         temperature = self.temperature
    #     if max_new_tokens is None:
    #         max_new_tokens = self.max_new_tokens
        
    #     if '<image>' in prompt:
    #         prompt = prompt.replace('<image>', '') # Not used for non vision models, this assumes that this class is always used for text models (as the vision model used is LLaVA and is implemented in a different class)
        
    #     messages = utils.get_messages(prompt)
    #     if additional_prompt is not None and '<ts>' in prompt:
    #         # import ipdb; ipdb.set_trace()
    #         assert type(additional_prompt) is np.ndarray, "additional_prompt must be a numpy ndarray"
    #         prompt = f"<|im_start|>system You are a helpful assistant.<|im_end|><|im_start|>user{prompt}<|im_end|><|im_start|>assistant"
    #         # build inputs via TSProcessor, then move each tensor to the target device
    #         time_series = [additional_prompt[:,t] for t in range(additional_prompt.shape[-1])]
    #         # import ipdb; ipdb.set_trace()
    #         inputs = self.ts_processor(text=[prompt], timeseries=time_series, padding='longest', return_tensors="pt").to(self.device)
    #     else:
    #         try:
    #             inputs = self.tokenizer.apply_chat_template(messages, add_generation_prompt=True, return_dict=True, return_tensors="pt").to(self.device)
    #         except:
    #             inputs = self.tokenizer(prompt, return_tensors="pt").to(self.device)
    #     # import ipdb; ipdb.set_trace()
    #     outputs = self.model.generate(**inputs, do_sample=True, temperature=temperature, top_k=self.top_k, top_p=self.top_p, num_beams=self.num_beams, 
    #                                 num_return_sequences=self.num_return_sequences, max_new_tokens=max_new_tokens, min_new_tokens=self.min_new_tokens, pad_token_id=self.tokenizer.eos_token_id)
    #     try:
    #         outputs = outputs[0][len(inputs[0]):] if not return_prompt else outputs[0]
    #     except:
    #         outputs = outputs[0][len(inputs['input_ids'][0]):] if not return_prompt else outputs[0]
    #     decoded_output = self.tokenizer.decode(outputs, skip_special_tokens=True)
        
    #     # Remove llama special words
    #     decoded_output = decoded_output.replace("assistant", "").replace("user", "").replace("system", "")

    #     return decoded_output
    
    def evaluate(self, ind):
        # transform GEPPY individual to Sympy expression
        raw_expr = gep.simplify(ind)
        sym_expr = parse_expr(str(raw_expr), transformations=standard_transformations)
        # print(sym_expr)
        # sym_expr now holds the Sympy form of the individual
        func = self.toolbox.compile(ind)
        n_correct = sum(func(*inp) == out for inp, out in zip(self.Input_data, self.Out_data))
        vars_used = {v.name for gene in ind for v in gene.kexpression}
        n_regs = len(vars_used - {"and_","or_","not_"})
        score = n_correct

        # print(score)

        # import ipdb; ipdb.set_trace()
        if self.use_llm:
            scoring_prompt = self.prompt.replace("{candidate_exprs}", str(sym_expr))
            llm_response = self.llm.generate(scoring_prompt)
            try:
                llm_score = float(re.findall(r"[-+]?\d*\.\d+|\d+", llm_response)[0])
            except (ValueError, IndexError):
                llm_score = 0.0
            try:
                complexity = sym_expr.count_ops()
            except:
                complexity = 0
            # import ipdb; ipdb.set_trace()
            score = n_correct + self.llm_weight*llm_score - 0*complexity

        return (score, n_regs)

    def infer_rules_for(self, target_node, nodes, transitions):
        # Build binary input/output lists
        self.Input_data = [[x_t[var] for var in nodes] for x_t, _ in transitions]
        self.Out_data = [x_tp1[target_node] for _, x_tp1 in transitions]
        # Build GEP primitive set
        pset = gep.PrimitiveSet("Main", input_names=nodes)
        pset.add_function(operator.and_, 2)
        pset.add_function(operator.or_, 2)
        pset.add_function(operator.not_, 1)
        # DEAP creator
        creator.create("FitnessMax", base.Fitness, weights=(1.0, -1.0))
        creator.create("Individual", gep.Chromosome, fitness=creator.FitnessMax)
        # Toolbox
        self.toolbox = gep.Toolbox()
        h = 5
        self.toolbox.register("gene_gen", gep.Gene, pset=pset, head_length=h)
        self.toolbox.register("individual", creator.Individual, gene_gen=self.toolbox.gene_gen, n_genes=1, linker=None)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("compile", gep.compile_, pset=pset)

        self.toolbox.register("evaluate", self.evaluate)
        self.toolbox.register("select", tools.selTournament, tournsize=3)
        # GEP operators
        self.toolbox.register("mut_uniform", gep.mutate_uniform, pset=pset, ind_pb=2/(2*h+1)); self.toolbox.pbs["mut_uniform"] = 0.5
        self.toolbox.register("mut_invert", gep.invert, pb=0.5)
        self.toolbox.register("mut_is_ts", gep.is_transpose, pb=0.5)
        self.toolbox.register("mut_ris_ts", gep.ris_transpose, pb=0.5)
        self.toolbox.register("mut_gene_ts", gep.gene_transpose, pb=0.5)
        self.toolbox.register("cx_1p", gep.crossover_one_point, pb=0.5)
        self.toolbox.register("cx_2p", gep.crossover_two_point, pb=0.5)
        self.toolbox.register("cx_gene", gep.crossover_gene, pb=0.5)
        # Run GEP
        # import ipdb; ipdb.set_trace()
        pop = self.toolbox.population(n=self.POP_SIZE)
        hof = tools.HallOfFame(self.hall_of_fame)
        pop2, log = gep_simple(pop, self.toolbox, n_generations=self.NUM_GENERATIONS, n_elites=self.hall_of_fame, hall_of_fame=hof, verbose=True)
        best = sorted(pop2, key=lambda ind: (-ind.fitness.values[0], ind.fitness.values[1]))[0]
        # import ipdb; ipdb.set_trace()
        expr = gep.simplify(best)
        return expr

    def fit(self, transitions, nodes, all_nodes):
        """Infer Boolean rules for each node given transitions and node order."""
        self.nodes = nodes
        self.all_nodes = all_nodes
        
        for node in nodes:
            self.inferred[node] = self.infer_rules_for(node, all_nodes, transitions)

        return self.inferred

    def predict(self, initial_state, steps=1):
        """Simulate Boolean network for a number of steps."""
        state = initial_state.copy()
        trajectory = [state.copy()]
        for _ in range(steps):
            new_state = {}
            for node, expr in self.inferred.items():
                f = lambdify(self.nodes, expr, 'numpy')
                args = [int(state[n]) for n in self.nodes]
                new_state[node] = bool(f(*args))
            state = new_state
            trajectory.append(state.copy())
        return trajectory

    def final_evaluate(self, inferred, dataset):
        train_transitions = dataset.data['train_transitions']
        test_transitions = dataset.data['test_transitions']
        # Use all input variables for evaluation
        nodes = self.all_nodes
        
        # Compute prediction metrics
        id_metrics = self.test_one_input(train_transitions, self.nodes, inferred)
        ood_metrics = self.test_one_input(test_transitions, self.nodes, inferred)
        # Compute average complexity of inferred expressions
        # import ipdb; ipdb.set_trace()
        complexities = [expr.count_ops() for expr in self.inferred.values()] if self.inferred else []
        id_metrics['complexity'] = float(np.mean(complexities)) if complexities else 0.0
        ood_metrics['complexity'] = id_metrics['complexity']

        print(f"ID Metrics: {id_metrics}")
        print(f"OOD Metrics: {ood_metrics}")

        return {'id': id_metrics, 'ood': ood_metrics}



    def test_one_input(self, transitions, nodes, inferred):
        """
        Evaluate inferred boolean rules against transitions and compute metrics.
        """
        # import ipdb; ipdb.set_trace()
        TP = FP = FN = TN = 0
        for node, expr in inferred.items():
            sympy_nodes = [Symbol(item) for item in self.all_nodes]
            f = lambdify(sympy_nodes, expr)
            for x_t, x_tp1 in transitions:
                pred = bool(f(**x_t))
                actual = bool(x_tp1[node])
                if   pred and actual: TP += 1
                elif pred and not actual: FP += 1
                elif not pred and actual: FN += 1
                else: TN += 1
        precision = TP/(TP+FP) if TP+FP>0 else 0
        recall    = TP/(TP+FN) if TP+FN>0 else 0
        f1        = 2*precision*recall/(precision+recall) if precision+recall>0 else 0
        accuracy  = (TP+TN)/(TP+FP+FN+TN) if TP+FP+FN+TN>0 else 0
        bm        = -1 + recall + (TN/(TN+FP) if TN+FP>0 else 0)
        return {'precision':precision,'recall':recall,'f1':f1,'accuracy':accuracy,'BM':bm}







