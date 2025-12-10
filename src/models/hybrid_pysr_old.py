import os
import fcntl
from utils import utils
import numpy as np
import random
from utils.data_utils import ExpressionTree
import torch
import re
from pysr import PySRRegressor
from odeformer.envs.generators import NodeList
import sympy
from scipy.integrate import solve_ivp
from odeformer.baselines.pysr_wrapper import PySRWrapper
from transformers import AutoTokenizer, AutoModelForCausalLM
# Assuming TSProcessor and snapshot_download are available in the environment
# from some_module import TSProcessor, snapshot_download
from odeformer.utils import timeout, MyTimeoutError


FITNESS_FAIL_SCORE = 9999.0
EPSILON = 1e-8
EARLY_STOP_TOLERANCE = 0.001
STALE_GEN_LIMIT = 3
MUTATION_STD_DEV = 0.02


class Hybrid_pysr(object):
    def __init__(self, model_name, device, dtype, cache_dir=None, **kwargs):
        self.model_name = model_name
        self.device = device
        self.dtype = dtype
        token = os.environ.get("HF_TOKEN", None)
        # import ipdb; ipdb.set_trace()
        if 'ChatTS-14B' in model_name:
            repo = "bytedance-research/ChatTS-14B"
            os.makedirs(cache_dir, exist_ok=True)
            model_path = snapshot_download(repo, cache_dir=cache_dir)
            # import ipdb; ipdb.set_trace()
            self.model = AutoModelForCausalLM.from_pretrained(model_path, trust_remote_code=True, device_map=self.device, torch_dtype=torch.float16)
            self.tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True, use_fast=False)
            self.ts_processor = TSProcessor(self.tokenizer, patch_size=16)
        else:
            self.tokenizer = AutoTokenizer.from_pretrained(model_name, torch_dtype=dtype, cache_dir=cache_dir, token=token)
            self.model = AutoModelForCausalLM.from_pretrained(model_name, torch_dtype=dtype, cache_dir=cache_dir, token=token, device_map='auto')
            if "tokenizer_pad" in kwargs:
                self.tokenizer.pad_token = kwargs["tokenizer_pad"]
            if "tokenizer_padding_side" in kwargs:
                self.tokenizer.padding_side = kwargs["tokenizer_padding_side"]

        self.model.eval()

        self.temperature = kwargs.get("temperature", 1.0)
        self.top_k = kwargs.get("top_k", 50)
        self.top_p = kwargs.get("top_p", 0.9)
        self.num_beams = kwargs.get("num_beams", 1)
        self.num_return_sequences = kwargs.get("num_return_sequences", 1)
        self.max_new_tokens = kwargs.get("max_new_tokens", 256)
        self.min_new_tokens = kwargs.get("min_new_tokens", 0)
        self.simplicity_weight = kwargs.get("simplicity_weight", 0.05)


        self.pysr_model = PySRWrapper(
            niterations=5,
            population_size=15,
            elementwise_loss="HuberLoss()",
            binary_operators=["+", "-", "*", "/", "^"],
            unary_operators=["exp", "log","sin", "cos"],
            verbosity=1,
            tournament_selection_n=10,
            optimize_hyperparams=False, # To avoid grid search logic in the wrapper
        )

        self.solve_config = {
            "t_span": (0, 10),
            "method": "BDF",       # switch to BDF to avoid LSODA prints
            "rtol": 1e-6,
            "atol": 1e-6,
            "first_step": 1e-6,
            "t_eval": np.linspace(0, 10, 150),
            "dense_output": False,
            }

        self._loss_function = timeout(seconds=60)(self._loss_function)

    def generate(self, prompt, additional_prompt=None, return_prompt=False, temperature=None, max_new_tokens=None):
        if temperature is None:
            temperature = self.temperature
        if max_new_tokens is None:
            max_new_tokens = self.max_new_tokens
        
        if '<image>' in prompt:
            prompt = prompt.replace('<image>', '') # Not used for non vision models, this assumes that this class is always used for text models (as the vision model used is LLaVA and is implemented in a different class)
        
        messages = utils.get_messages(prompt)
        if additional_prompt is not None and '<ts>' in prompt:
            # import ipdb; ipdb.set_trace()
            assert type(additional_prompt) is np.ndarray, "additional_prompt must be a numpy ndarray"
            prompt = f"<|im_start|>system You are a helpful assistant.<|im_end|><|im_start|>user{prompt}<|im_end|><|im_start|>assistant"
            # build inputs via TSProcessor, then move each tensor to the target device
            time_series = [additional_prompt[:,t] for t in range(additional_prompt.shape[-1])]
            # import ipdb; ipdb.set_trace()
            inputs = self.ts_processor(text=[prompt], timeseries=time_series, padding='longest', return_tensors="pt").to(self.device)
        else:
            try:
                inputs = self.tokenizer.apply_chat_template(messages, add_generation_prompt=True, return_dict=True, return_tensors="pt").to(self.device)
            except:
                inputs = self.tokenizer(prompt, return_tensors="pt").to(self.device)
        # import ipdb; ipdb.set_trace()
        outputs = self.model.generate(**inputs, do_sample=True, temperature=temperature, top_k=self.top_k, top_p=self.top_p, num_beams=self.num_beams, 
                                    num_return_sequences=self.num_return_sequences, max_new_tokens=max_new_tokens, min_new_tokens=self.min_new_tokens, pad_token_id=self.tokenizer.eos_token_id)
        try:
            outputs = outputs[0][len(inputs[0]):] if not return_prompt else outputs[0]
        except:
            outputs = outputs[0][len(inputs['input_ids'][0]):] if not return_prompt else outputs[0]
        decoded_output = self.tokenizer.decode(outputs, skip_special_tokens=True)
        
        # Remove llama special words
        decoded_output = decoded_output.replace("assistant", "").replace("user", "").replace("system", "")

        return decoded_output


    
    def _loss_function(self, weights):
        """
        Compute weighted fitness using symbolic regression + LLM + complexity.
        """
        if len(weights) != 3 or not np.isclose(weights.sum(), 1.0, atol=EPSILON):
            return FITNESS_FAIL_SCORE

        w1, w2, w3 = weights
        try:
            # fit the symbolic regression model
            # import ipdb; ipdb.set_trace()
            self.pysr_model.fit(self.times, self.trajectories, verbose=False, sort_candidates=False)
            best = self.pysr_model.get_best()
            eq_str = [item.equation for item in best]
            eq_str
            best_list = [item.sympy_format for item in best]
            if not eq_str:
                return FITNESS_FAIL_SCORE

            var_symbols = [sympy.symbols(f"x_{i}") for i in range(self.trajectories.shape[1])]
            # import ipdb; ipdb.set_trace()
            callable_fn = lambda t, x: np.array([f(*x) for f in [sympy.lambdify(var_symbols, eq, 'numpy') for eq in best_list]])
            init_conds = self.trajectories[0]
            # suppress invalid-value-in-power warnings during integration
            with np.errstate(invalid='ignore', divide='ignore'):
                # import ipdb; ipdb.set_trace()
                sol = solve_ivp(callable_fn, **self.solve_config, y0=init_conds)
            
            y_pred = sol.y.T
            mse = float(np.mean((self.trajectories - y_pred) ** 2))

            # import ipdb; ipdb.set_trace()
            # LLM evaluation via generate
            scoring_prompt = self.prompt.replace("{candidate_exprs}", str(eq_str))
            llm_response = self.generate(scoring_prompt)
            try:
                llm_score = float(re.findall(r"[-+]?\d*\.\d+|\d+", llm_response)[0])
            except (ValueError, IndexError):
                llm_score = 0.0

            # complexity
            complexity = sum([item.count_ops() for item in best_list])

            return w1 * mse + w2 * (1 - llm_score/10) + w3 * complexity
        except Exception:
            return FITNESS_FAIL_SCORE

    def fit(self, trajectories, times=None, prompt_text="", generations=50, pop_size=15, mutation_rate=0.3):
        # store inputs
        self.trajectories = trajectories
        self.times = times
        self.prompt_text = prompt_text

        # initialize population of weight vectors
        pop = [np.random.dirichlet([3, 2, 1]) for _ in range(pop_size)]
        stale_gens = 0
        prev_best = float('inf')

        for gen in range(generations):
            fitness_vals = np.array([self._loss_function(ind) for ind in pop])
            order = fitness_vals.argsort()
            best_fit = float(fitness_vals[order[0]])
            avg_fit = float(fitness_vals.mean())

            # early stopping check
            improvement = (prev_best - best_fit) / (abs(prev_best) + EPSILON)
            if improvement < EARLY_STOP_TOLERANCE:
                stale_gens += 1
            else:
                stale_gens = 0
            prev_best = best_fit
            if stale_gens >= STALE_GEN_LIMIT:
                print(f"[GA] Early stopping at generation {gen} – fitness plateaued.")
                break
            
            # import ipdb; ipdb.set_trace()
            # selection and reproduction
            elites = [pop[i] for i in order[:pop_size // 2]]
            offspring = []
            for _ in range(pop_size // 2):
                p1, p2 = random.sample(elites, 2)
                child = (p1 + p2) / 2
                if random.random() < mutation_rate:
                    child += np.random.normal(0, MUTATION_STD_DEV, size=3)
                child = np.clip(child, 0, 1)
                child = child / child.sum() if child.sum() > 0 else np.array([1, 0, 0])
                offspring.append(child)
            pop = elites + offspring

        # final selection of best weight vector
        fitness_vals = np.array([self._loss_function(ind) for ind in pop])
        best_idx = int(fitness_vals.argmin())
        return pop[best_idx], float(fitness_vals[best_idx])



