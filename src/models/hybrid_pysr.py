import os
from utils import utils
import numpy as np
import random
from utils.data_utils import ExpressionTree

from gplearn.genetic import SymbolicRegressor
import torch
import re
from odeformer.envs.generators import NodeList
import sympy
from scipy.integrate import solve_ivp
from sympy import sympify, lambdify
from transformers import AutoTokenizer, AutoModelForCausalLM
# Assuming TSProcessor and snapshot_download are available in the environment
# from some_module import TSProcessor, snapshot_download
from odeformer.utils import timeout, MyTimeoutError
import pickle
import tempfile


class Hybrid_DE(object):
    def __init__(self, model_name, device, dtype, cache_dir=None, **kwargs):
        self.model_name = model_name
        self.device = device
        self.dtype = dtype
        token = os.environ.get("HF_TOKEN", None)
        # initialize LLM model and tokenizer as before
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
        
        # Parsimony coefficient for gplearn
        self.complexity_weight = kwargs.get("complexity_weight", 0.1)

        # Store data for use in loss function
        self.current_times = None
        self.current_trajectories = None
        self.current_prompt = None

        # Initialize gplearn SymbolicRegressor for each variable
        gp_kwargs = {
            "population_size": kwargs.get("population_size", 100),
            "generations": kwargs.get("niterations", 20),
            "function_set": ['add', 'sub', 'mul', 'div', 'sin', 'cos', 'exp', 'log'],
            "parsimony_coefficient": self.complexity_weight,
            "metric": 'mean squared error',
            "p_crossover": 0.9,
            "p_subtree_mutation": 0.01,
            "p_hoist_mutation": 0.01,
            "p_point_mutation": 0.01,
            "max_samples": 1.0,
            "verbose": 1,
            "random_state": 0,
        }
        self.gp_kwargs = gp_kwargs
        self.gp_models = []  # one SymbolicRegressor per variable

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
            

    def fit(self, trajectories, times=None, prompt_text="", **kwargs):
        """
        Fit the hybrid PySR model to trajectory data.
        
        Args:
            trajectories: numpy array of shape (n_timepoints, n_variables)
            times: numpy array of time points
            prompt_text: prompt template for LLM evaluation (should contain {candidate_exprs})
            **kwargs: additional arguments passed to PySR
        """
        # Compute time points
        self.current_trajectories = trajectories
        self.current_times = times if times is not None else np.linspace(0, 10, trajectories.shape[0])
        # Estimate derivatives via finite differences
        dt = np.gradient(self.current_times)
        derivatives = np.gradient(trajectories, axis=0) / dt[:, None]
        # Prepare training inputs
        X = trajectories
        # Fit one SymbolicRegressor per variable using stored gp_kwargs
        from sklearn.utils.validation import check_X_y
        for i in range(trajectories.shape[1]):
            model = SymbolicRegressor(**self.gp_kwargs)
            # Monkey-patch validate method for sklearn compatibility
            if not hasattr(model, '_validate_data'):
                def _validate_data(self, X, y=None, **kwargs):
                    return check_X_y(X, y, **kwargs)
                model._validate_data = _validate_data.__get__(model, model.__class__)
            model.fit(X, derivatives[:, i])
            self.gp_models.append(model)
        return self.get_best_equations()

    def get_best_equations(self):
        # Return expressions for each variable
        equations = {i: [model._program.__str__()] for i, model in enumerate(self.gp_models)}
        return equations

    def predict(self, times, initial_conditions=None):
        """Make predictions using the best found equations."""
        # Build ODE system from learned expressions
        # Prepare lambdified functions
        funcs = []
        for model in self.gp_models:
            expr = sympify(model._program.__str__())
            f = lambdify([*['x'+str(i) for i in range(len(self.gp_models))]], expr, 'numpy')
            funcs.append(f)
        # Define ODE rhs
        def rhs(t, y):
            return [func(*y) for func in funcs]
        # Set initial conditions
        y0 = initial_conditions if initial_conditions is not None else self.current_trajectories[0]
        sol = solve_ivp(rhs, (times[0], times[-1]), y0, t_eval=times)
        return sol.y.T





