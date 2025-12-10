import os
import sys
import json
import copy
import time
import datetime
import warnings
warnings.filterwarnings("ignore", message="invalid value encountered in scalar power")
import signal
import cProfile
from typing import Dict, Tuple, List, Any
from collections.abc import Callable

# Import juliacall BEFORE torch to avoid segfault
# Set environment variable to handle signals properly
os.environ['PYTHON_JULIACALL_HANDLE_SIGNALS'] = 'yes'
try:
    import juliacall
except ImportError:
    print("Warning: juliacall not available, Julia features will be disabled")
    juliacall = None

import hydra
import torch
import matplotlib.pyplot as plt
import numpy as np
from transformers import set_seed
from omegaconf import OmegaConf, DictConfig, listconfig
from sklearn.metrics import r2_score
import argparse

from utils import utils
from verifiers import get_verifier
from datasets import get_dataset
from models import load_model

from mloggers import ConsoleLogger, FileLogger, MultiLogger, LogLevel

import ipdb


class Pipeline(object):
    def __init__(self, cfg: DictConfig) -> None:
        # ipdb.set_trace()
        self.cfg = cfg
        # -------- initializing settings --------
        # Output setup
        self.root_dir = cfg.get("root", os.getcwd())
        self.output_dir = cfg.get("output_dir", "output")
        model_folder_name = cfg.model.name.strip()
        if "/" in model_folder_name:
            model_folder_name = model_folder_name.split("/")[-1]
        experiment_folder_name = os.path.join(cfg.model.name.split("/")[-1], cfg.experiment.input_type ,cfg.experiment.symbolic_expression.name)

        expr_idx = cfg.experiment.symbolic_expression.get("idx", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        expr_idx = "Expr_"+str(expr_idx)

        self.output_path = os.path.join(self.root_dir, self.output_dir, experiment_folder_name, expr_idx + "/")
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        # Logger setup
        cfg.logger.run_id = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        loggers_list = cfg.logger.loggers
        log_level = LogLevel[cfg.logger.get("level", "INFO")]
        loggers = []
        for logger in loggers_list:
            if logger == "console":
                loggers.append(ConsoleLogger(default_priority=log_level))
            elif logger == "file":
                loggers.append(FileLogger(os.path.join(self.output_path, 'log.json'), default_priority=log_level))
            elif logger == "":
                pass
            else:
                print(f'[WARNING] Logger "{logger}" is not supported')
        self.logger = MultiLogger(loggers, default_priority=log_level)
        self.logger.info(f"Project root: {self.root_dir}.")
        self.logger.info(f"Logging to {self.output_path}.")
        job_id = utils.get_job_id()
        self.logger.info(f"Slurm job ID: {job_id}.") if job_id is not None else None

        # Redirect warnings to logger
        warnings.filterwarnings("default")
        warnings.showwarning = lambda *args, **kwargs: self.logger.warning(str(args[0]))

        self.interval_save = cfg.experiment.get("interval_save", 1)
        
        # RNG setup
        if not hasattr(cfg, "seed") or cfg.seed is None or cfg.seed == -1:
            self.cfg.seed = np.random.randint(0, np.iinfo(np.int32).max)
            self.logger.info(f"Seed not specified, using random seed: {self.cfg.seed}.")
        else:
            self.logger.info(f"Using seed: {self.cfg.seed}.")

        np.random.seed(self.cfg.seed)
        torch.manual_seed(self.cfg.seed)
        torch.cuda.manual_seed_all(self.cfg.seed) if torch.cuda.is_available() else None
        set_seed(self.cfg.seed)

        self.max_retries = self.cfg.experiment.symbolic_expression.get("max_retries", 10)

        if torch.cuda.is_available():
            torch.cuda.init()

        if cfg.device == "auto":
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = torch.device(cfg.device)

        if cfg.get("use_bfloat16", False):
            self.dtype = torch.bfloat16 if torch.cuda.is_available() and torch.cuda.is_bf16_supported() else torch.float16
        else:
            self.dtype = torch.float16

        self.logger.info(f"Using device: {self.device} with dtype: {self.dtype}.")
        if torch.cuda.is_available() and ('cuda' in cfg.device or 'auto' in cfg.device):
            self.logger.info(f"Device name: {torch.cuda.get_device_name()}.")
        
        self.cache_dir = self.cfg.model.get("cache_dir", os.environ.get("HF_HOME", None))
        if self.cache_dir == "":
            self.cache_dir = os.environ.get("HF_HOME", None)
        
        if self.cache_dir is not None:
            os.environ['HF_HOME'] = self.cache_dir 
            os.environ['TRANSFORMERS_CACHE'] = os.environ['HF_HOME']
        self.logger.info(f"Cache dir: {os.environ.get('HF_HOME', None)}.")
        
        # build dataset
        self.dataset = get_dataset(cfg, logger=self.logger)

        self.verifier = get_verifier(cfg, self.dataset, self.logger, self.device)

        self.additional_prompt = self.verifier.additional_prompt
        # Model settings
        self.build_model(cfg)

        self.prompt_path = cfg.experiment.get("textual_input", None)
        with open(self.prompt_path, 'r') as f:
            self.prompt = json.load(f)
        
        # import ipdb; ipdb.set_trace()

        self.model.prompt = self.verifier.prepare_prompt(self.prompt['judge_prompt'])
        self.model.additional_prompt = self.verifier.additional_prompt

        # import ipdb; ipdb.set_trace()


    
    def build_model(self, cfg) -> None:
        # import ipdb; ipdb.set_trace()
        self.model_name = cfg.model.name
        self.model = None

        self.tokenizer_pad = cfg.model.tokenizer_pad
        self.tokenizer_padding_side = cfg.model.tokenizer_padding_side

        self.max_new_tokens = cfg.model.max_new_tokens
        self.top_p = cfg.model.top_p
        self.top_k = cfg.model.top_k
        self.num_beams = cfg.model.num_beams

        self.generation_tokens = self.cfg.experiment.max_output_tokens if hasattr(self.cfg.experiment, "max_output_tokens") else 512
        self.temperature = cfg.model.temperature
        if cfg.model.temperature_schedule:
            self.temperature_scheduler = torch.optim.lr_scheduler.ExponentialLR(torch.optim.Adam([torch.tensor(self.temperature)], lr=1), gamma=cfg.model.temperature_schedule_gamma)
        else:
            self.temperature_scheduler = None
            
        model_args = {
            "temperature": self.temperature,
            "top_p": self.top_p,
            "top_k": self.top_k,
            "num_beams": self.num_beams,
            "max_length": self.max_new_tokens,
            "min_length": 0,
            "tokenizer_pad": self.tokenizer_pad,
            "tokenizer_padding_side": self.tokenizer_padding_side,
            "seed": self.cfg.seed,
            "api_key_path": os.path.join(self.root_dir, cfg.model.api_key_path) if hasattr(cfg.model, "api_key_path") else None,
            "organization_id_path": os.path.join(self.root_dir, cfg.model.organization_id_path) if hasattr(cfg.model, "organization_id_path") else None,
        }
        if torch.cuda.is_available() and 'A100' in torch.cuda.get_device_name(0):
            model_args['attn_implementation'] = 'flash_attention_2'
            model_args['use_flash_attn'] = True
            self.logger.info("Using Flash Attention 2")

        self.model = load_model(self.model_name, self.device, self.dtype, self.cache_dir, model_args, exp_name=cfg.experiment.symbolic_expression.name)
        self.logger.info(f"Loaded Model: {cfg.model.name}, visual_inputs: {cfg.model.visual}, time_series_embedding: {cfg.model.time_series}.")
        self.model.use_llm = cfg.use_llm
        self.model.llm_weight = cfg.llm_weight
        print(f"Using LLM: {cfg.use_llm}, weight: {cfg.llm_weight}.")
        # import ipdb; ipdb.set_trace()



    def run(self) -> None:
        if "DE" in self.cfg.experiment.symbolic_expression.name:
            times = self.verifier.solve_config['t_eval']
            fitted_results = self.model.fit(self.verifier.time_series, times)
            metrics = self.model.final_evaluate(fitted_results, self.dataset)
        elif "BN" in self.cfg.experiment.symbolic_expression.name:
            all_nodes = self.verifier.dataset.var_list
            nodes = self.verifier.dataset.trans_vars
            fitted_results = self.model.fit(self.verifier.time_series, nodes, all_nodes)
            metrics = self.model.final_evaluate(fitted_results, self.dataset)
        print(f"Fitted Results: {fitted_results}")
        print(f"Metrics: {metrics}")
        # import ipdb; ipdb.set_trace()
        np.save(self.output_path + f"final_results.npy""", {"exprs": fitted_results, "metrics": metrics})
        self.logger.info(f"Final Results Saved.")
        


@hydra.main(version_base=None, config_path="conf", config_name="config")
def main(cfg: DictConfig) -> None:
    # Set environment variable to handle signals properly and import juliacall before torch
    os.environ['PYTHON_JULIACALL_HANDLE_SIGNALS'] = 'yes'
    try:
        import juliacall
    except ImportError:
        print("Warning: juliacall not available, Julia features will be disabled")
        juliacall = None
        
    workspace = Pipeline(cfg)
    workspace.run()

def dump_profile():
    profiler.disable()
    job_id = utils.get_job_id()
    print(f"Dumping profile to {os.path.join(os.getcwd(), 'profiles', 'profile')}_{job_id if job_id is not None else 'local'}")
    if not os.path.exists("./profiles"):
            os.makedirs("./profiles")
    profiler.dump_stats(f"./profiles/profile_{job_id if job_id is not None else 'local'}")

def signal_handler(sig, frame):
    # Ignore warnings, as otherwise we break the logger
    warnings.filterwarnings("ignore")
    dump_profile()
    print(f"Detecting signal {sig}. Dumping profile to {os.path.join(os.getcwd(), 'profiles', 'profile')}_{job_id if job_id is not None else 'local'}")
    sys.stdout.flush()
    if sig == signal.SIGTERM or sig == signal.SIGINT:
        sys.exit(1)

if __name__ == "__main__":

    # Run full profiler if env variable PROFILE is set
    do_profile = os.environ.get("PROFILE", False)
    print("Initializing profiler.")
    print("Profile will only be created if the code fails or is terminated.") if not do_profile else print("Profile will be created.")
    job_id = utils.get_job_id()

    # Set termination signal handlers to dump profile when terminated by SLURM
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGCONT, signal_handler)
    
    # Setup profiler
    global profiler
    profiler = cProfile.Profile()
    profiler.enable()

    try:
        main()
    except Exception as e:
        # Catch exceptions and dump profile
        print("Caught exception in main.")
        print(e)
        dump_profile()
        sys.exit(2)

    if do_profile:
        dump_profile()
    print("Finished main.")