import os
import sys
import json
import copy
import time
import datetime
import warnings
warnings.filterwarnings("ignore")
import logging
logging.getLogger().setLevel(logging.ERROR)
import signal
import cProfile
from typing import Dict, Tuple, List, Any
from collections.abc import Callable
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

# from mloggers import ConsoleLogger, FileLogger, MultiLogger, LogLevel
# Create simple logger classes to avoid import errors and warnings

class LogLevel:
    DEBUG = 0
    INFO = 1
    WARNING = 2
    ERROR = 3

class ConsoleLogger:
    def __init__(self, default_priority=LogLevel.INFO):
        self.default_priority = default_priority
    
    def info(self, msg):
        if self.default_priority <= LogLevel.INFO:
            print(f"[INFO] {msg}")
    
    def warning(self, msg):
        if self.default_priority <= LogLevel.WARNING:
            print(f"[WARNING] {msg}")
    
    def error(self, msg):
        if self.default_priority <= LogLevel.ERROR:
            print(f"[ERROR] {msg}")

class FileLogger:
    def __init__(self, filepath, default_priority=LogLevel.INFO):
        self.filepath = filepath
        self.default_priority = default_priority
    
    def info(self, msg):
        pass  # Suppress file logging to avoid warnings
    
    def warning(self, msg):
        pass  # Suppress file logging to avoid warnings
    
    def error(self, msg):
        pass  # Suppress file logging to avoid warnings

class MultiLogger:
    def __init__(self, loggers, default_priority=LogLevel.INFO):
        self.loggers = loggers
        self.default_priority = default_priority
    
    def info(self, msg):
        for logger in self.loggers:
            logger.info(msg)
    
    def warning(self, msg):
        for logger in self.loggers:
            logger.warning(msg)
    
    def error(self, msg):
        for logger in self.loggers:
            logger.error(msg)

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
        log_level = getattr(LogLevel, cfg.logger.get("level", "INFO"), LogLevel.INFO)
        loggers = []
        for logger in loggers_list:
            if logger == "console":
                loggers.append(ConsoleLogger(default_priority=log_level))
            elif logger == "file":
                # Use our custom FileLogger that suppresses warnings
                loggers.append(FileLogger(os.path.join(self.output_path, 'log.json'), default_priority=log_level))
            elif logger == "":
                pass
            else:
                # Suppress this warning message too
                pass
        self.logger = MultiLogger(loggers, default_priority=log_level)
        self.logger.info(f"Project root: {self.root_dir}.")
        self.logger.info(f"Logging to {self.output_path}.")
        job_id = utils.get_job_id()
        self.logger.info(f"Slurm job ID: {job_id}.") if job_id is not None else None

        # Redirect warnings to logger - keep warnings suppressed
        warnings.filterwarnings("ignore")
        # Suppress logging warnings from libraries
        logging.getLogger("gplearn").setLevel(logging.ERROR)
        logging.getLogger("ensembling").setLevel(logging.ERROR)
        logging.getLogger("numpy").setLevel(logging.ERROR)
        logging.getLogger().setLevel(logging.ERROR)

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
        # build base prompts
        self.build_prompts(cfg)
        # Model settings
        self.build_model(cfg)
        # Verifier
        # import ipadb; ipdb.set_trace()
        self.verifier = get_verifier(cfg, self.dataset, self.logger, self.device)

        self.hybrid = cfg.hybrid
        self.use_llm = cfg.use_llm
        self.llm_weight = cfg.llm_weight

    def build_prompts(self, cfg) -> None:

        if cfg.experiment.input_type == "textual":
            with open(os.path.join(self.root_dir, cfg.experiment.textual_input), "r") as f:
                self.prompts = json.load(f)
        elif cfg.experiment.input_type == "visual":
            with open(os.path.join(self.root_dir, cfg.experiment.visual_input), "r") as f:
                self.prompts = json.load(f)
        elif cfg.experiment.input_type == "time_series":
            with open(os.path.join(self.root_dir, cfg.experiment.time_series_input), "r") as f:
                self.prompts = json.load(f)
        else:
            self.logger.error(f"Input type {cfg.expriment.input_type} not supported.")
            exit(1)

        self.seed_prompt = self.prompts['seed_prompt']
        self.base_prompt = self.prompts['base_prompt']
    
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

        self.iterations = cfg.experiment.symbolic_expression.get("iterations", 10)


        self.generation_tokens = self.max_new_tokens
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

        # import ipdb; ipdb.set_trace()
        self.model = load_model(self.model_name, self.device, self.dtype, self.cache_dir, model_args)
        self.logger.info(f"Loaded Model: {cfg.model.name}, visual_inputs: {cfg.model.visual}, time_series_embedding: {cfg.model.time_series}.")
    
    def generate(self, prompt, additional_prompt=None, max_tries=10, early_stop=False):
        
        # max_tries = self.cfg.experiment.max_retri
        functions = []
        reasonings = []

        start_time = time.perf_counter()
        with torch.inference_mode():
            # import ipdb; ipdb.set_trace()
            for i in range(max_tries):
                # Generate seed functions using the model
                self.logger.info(f"Attempt {i+1} of {max_tries} to generate functions.")
                
                model_output = self.model.generate(prompt, additional_prompt=additional_prompt, temperature=self.temperature, max_new_tokens=self.generation_tokens)
                # import ipdb; ipdb.set_trace()
                # clean_prompt = f"Check if the outputs are valid functions. Make sure the variables are chosen from the list {str(self.verifier.var_list)}. Refine the current output: {model_output}. Only return the valid functions. Do not add any additional text."
                # self.model.generate(clean_prompt, additional_prompt=additional_prompt, temperature=self.temperature, max_new_tokens=self.generation_tokens)
                # import ipdb; ipdb.set_trace()
                # self.logger.info(f"Model output: {model_output}.")
                # import ipdb; ipdb.set_trace()
                new_functions, new_reasons = self.verifier.parse_model_output(model_output)
                self.logger.info("Parsed model outputs: " + str(new_functions))


                # import ipdb; ipdb.set_trace()
                functions += new_functions
                reasonings += new_reasons

                if len(functions) > 50:
                    break
                # if (early_stop and len(functions) > 1) or len(functions) > 25:
                #     break

                # import ipdb; ipdb.set_trace()
                # valid_functions = self.verifier.check_validity(new_functions)
                # import ipdb; ipdb.set_trace()
                # if len(valid_functions) != 0:
                #     functions.update(new_functions)

                # if len(functions) > 5:
                #     break
        end_time = time.perf_counter()
        # import ipdb; ipdb.set_trace()

        self.logger.info(f"Generated {len(functions)} functions: {functions}.")
        # import ipdb; ipdb.set_trace()

        return list(zip(functions, reasonings)), end_time - start_time


    def inialization(self, ):
        # test = 'c*x_0 - c*x_1*x_2 | c*x_1 + c*x_0*x_2 | c*x_2 + c*x_0 * x_1'
        # exprs = self.verifier.dataset.str_to_sympy(test)
        # fit_exprs = self.verifier.fit(exprs, self.verifier.init_conds, self.verifier.time_series)
        # score = utils.r2_score(self.verifier.time_series, self.verifier.generate([fit_exprs], [self.verifier.init_conds])[0])
        # import ipdb; ipdb.set_trace()
        # Generate seed functions
        self.seed_functions = {}
        # min_seed_functions = self.cfg.experiment.min_seed_functions
        min_seed_functions = 1
        gen_time = 0

        # prepare seed prompt
        self.seed_prompt = self.verifier.prepare_prompt(self.seed_prompt)
        print("seed_prompt:", self.seed_prompt)
        # import ipdb; ipdb.set_trace()
        self.seed_functions, gen_time = self.generate(self.seed_prompt, additional_prompt=self.verifier.additional_prompt, max_tries=5, early_stop=False)
        
        # import ipdb; ipdb.set_trace()
        assert len(tuple(self.seed_functions)) >= min_seed_functions, f"Could not generate {min_seed_functions} seed functions. Generated {len(tuple(self.seed_functions))} seed functions."
        # import ipdb; ipdb.set_trace()
        
        self.verifier.update_candidates(self.seed_functions, hybrid=self.hybrid)

        # import ipdb; ipdb.set_trace()

        self.logger.info(f"Current functions: {self.verifier.candidate_expressions}.")
        self.logger.info(f"Succesfully generated {len(tuple(self.seed_functions))} seed functions in {gen_time} seconds.")

        # import ipdb; ipdb.set_trace()
        # Results json
        self.results = {
            "experiment_name": self.cfg.experiment.symbolic_expression.name,
            "seed": self.cfg.seed,
            "train_points": self.verifier.time_series,
            "test_points": self.verifier.time_series,
            "candidates": self.verifier.candidate_expressions,
            "iterations": 0,
            "tries_per_iteration": [],
            "epoch_scores": {},
            "times": {
                "iteration": [],
                "seed_function_generation": gen_time,
                "generation_per_iteration": [],
                "optimization_per_iteration": [],
            }
        }

        # Save config
        with open(self.output_path + "config.yaml", "w") as f:
            OmegaConf.save(self.cfg, f)
    

    def iterative_refinement(self, iterations=None):
        # import ipdb; ipdb.set_trace()
        main_timer_start = time.perf_counter()
        # Start the main loop
        for i in range(iterations):
            start_time = time.perf_counter()
            self.logger.info(f"Round {i+1}.")

            new_prompt = self.verifier.prepare_prompt(self.base_prompt)
            print("new prompt:", new_prompt)

            iter_functions, gen_time = self.generate(new_prompt, additional_prompt=self.verifier.additional_prompt, max_tries=5, early_stop=False)
            
            self.verifier.update_candidates(iter_functions, hybrid=self.hybrid)

            # import ipdb; ipdb.set_trace()

            self.logger.info(f"Current functions: {self.verifier.candidate_expressions.head()}.")
            self.logger.info(f"Succesfully generated {len(tuple(iter_functions))} functions in {gen_time} seconds.")

            # Handle temperature schedule
            if self.temperature_scheduler is not None:
                self.temperature = self.temperature_scheduler.get_last_lr()[0]
                self.logger.info(f"Temperature: {self.temperature}.")
                self.results["temperatures"].append(self.temperature)
                self.temperature_scheduler.step()
            try:
                self.results['epoch_scores'][i] = self.verifier.candidate_expressions.iloc[0]
            except:
                pass
            
            # import ipdb; ipdb.set_trace()
            # save checkpoint
            if i % self.interval_save == 0:
                self.logger.info(f"Checkpoint {i}. Saving results.")
                checkpoint_timer_end = time.perf_counter()
                self.results["times"]["total"] = checkpoint_timer_end - main_timer_start
                self.results["candidates"] = self.verifier.candidate_expressions
                # import ipdb; ipdb.set_trace()
                try:
                    np.save(self.output_path + f"results_checkpoint_{i}.npy", self.results)
                except:
                    import dill
                    with open(self.output_path + f"results_checkpoint_{i}.pkl", "wb") as f:
                        dill.dump(self.results, f)
                    
                self.logger.info(f"Checkpoint {i} saved.")
            
            valid = self.verifier.check_tolerance()

            if valid or i == (iterations - 1):
                # expr, function, score = self.verifier.get_best_expr()
                self.results["last_found_at"] = i
                timer_end = time.perf_counter()
                self.results["times"]["total"] = timer_end - main_timer_start
                self.results["candidates"] = self.verifier.candidate_expressions
                # import ipdb; ipdb.set_trace()

                if i == iterations - 1:
                    self.logger.info(f"Iterations Finished.")

                return valid
            
        return False


    def run(self) -> None:
        # generate initial proposals
        # import ipdb; ipdb.set_trace()
        self.inialization()
        # import ipdb; ipdb.set_trace()
        # check if the inital proposals meet the tolerance
        # import ipdb; ipdb.set_trace()
        if self.verifier.check_tolerance():
            # expr, function, score = self.verifier.get_best_expr()
            self.logger.info(f"Initial proposals already meet the tolerance.")
            self.results["epoch_scores"][0] = self.verifier.candidate_expressions.iloc[0]
            # import ipdb; ipdb.set_trace()
            valid = True
        else:
            valid = self.iterative_refinement(iterations=self.iterations)
    
        
        if valid:
            self.logger.info(f"Found valid function.")
        else:
            self.logger.info(f"Did not find valid function.")

        # import ipdb; ipdb.set_trace()
        if valid and self.verifier.OOD:
            ood_results = self.verifier.check_OOD()
            self.results["OOD"] = ood_results

        # import ipdb; ipdb.set_trace()
        np.save(self.output_path + f"final_results.npy""", self.results)
        self.logger.info(f"Final Results Saved.")
        
            



@hydra.main(version_base=None, config_path="conf", config_name="config")
def main(cfg: DictConfig) -> None:
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
    logging.getLogger().setLevel(logging.ERROR)
    dump_profile()
    print(f"Detecting signal {sig}. Dumping profile to {os.path.join(os.getcwd(), 'profiles', 'profile')}_{job_id if job_id is not None else 'local'}")
    sys.stdout.flush()
    if sig == signal.SIGTERM or sig == signal.SIGINT:
        sys.exit(1)

if __name__ == "__main__":
    # Suppress all warnings and logging messages
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    
    # Suppress specific library warnings
    import os
    os.environ['PYTHONWARNINGS'] = 'ignore'
    
    # Suppress numpy warnings
    import numpy as np
    np.seterr(all='ignore')

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