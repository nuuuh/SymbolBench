import os
import fcntl
from utils import utils
import numpy as np
import torch

from transformers import AutoTokenizer, AutoModelForCausalLM
from huggingface_hub import snapshot_download

os.environ["HF_TOKEN"] = "your_HF_token_here"


class TSProcessor:
    def __init__(self, tokenizer, patch_size):
        self.tokenizer = tokenizer
        self.patch_size = patch_size
    def __call__(self, text, timeseries, padding=True, return_tensors="pt"):
        token_outputs = self.tokenizer(text, padding=padding, return_tensors=return_tensors)
        ts_list = []
        for ts in timeseries:
            # import ipdb; ipdb.set_trace()
            arr = np.array(ts, dtype=np.float16)

            # Replace inf values with 0 and create a mask for valid values
            valid_mask = np.isfinite(arr).astype(np.float16)  # Mask where finite values are 1, and inf/nan are 0
            arr = np.where(valid_mask, arr, 0)  # Replace inf/nan values with 0

            L = arr.shape[-1]
            pad_len = (-L) % self.patch_size
            if pad_len:
                arr_padded = np.pad(arr, (0, pad_len), mode='constant')
                valid_mask_padded = np.pad(valid_mask, (0, pad_len), mode='constant')
            else:
                arr_padded = arr
                valid_mask_padded = valid_mask

            mask = np.stack((valid_mask_padded, np.zeros_like(valid_mask_padded)), axis=1)  # shape: (seq_length, num_features)
            feature = np.stack((arr_padded, mask[:, 0]), axis=1)  # shape: (seq_length, num_features)
            ts_list.append(torch.tensor(feature, dtype=torch.float16))
        # import ipdb; ipdb.set_trace()
        ts_tensor = torch.stack(ts_list, dim=0)  # shape: (batch, sequence_length, num_features)
        ts_tensor = ts_tensor.to(torch.float16)  # cast to model dtype
        token_outputs['timeseries'] = ts_tensor
        return token_outputs



class HuggingFaceModel(object):
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
            # ChatTS-14B doesn't have safetensors, so we don't specify use_safetensors
            self.model = AutoModelForCausalLM.from_pretrained(model_path, trust_remote_code=True, device_map=self.device, torch_dtype=torch.float16)
            self.tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True, use_fast=False)
            self.ts_processor = TSProcessor(self.tokenizer, patch_size=16)
        else:
            self.tokenizer = AutoTokenizer.from_pretrained(model_name, torch_dtype=dtype, cache_dir=cache_dir, token=token)
            self.model = AutoModelForCausalLM.from_pretrained(model_name, torch_dtype=dtype, cache_dir=cache_dir, token=token, device_map=self.device)
            if "tokenizer_pad" in kwargs:
                self.tokenizer.pad_token = kwargs["tokenizer_pad"]
            if "tokenizer_padding_side" in kwargs:
                self.tokenizer.padding_side = kwargs["tokenizer_padding_side"]

        self.model.eval()
        # self.model = self.model.bfloat16()

        # import ipdb; ipdb.set_trace()

        self.temperature = kwargs.get("temperature", 1.0)
        self.top_k = kwargs.get("top_k", 50)
        self.top_p = kwargs.get("top_p", 0.9)
        self.num_beams = kwargs.get("num_beams", 1)
        self.num_return_sequences = kwargs.get("num_return_sequences", 1)
        self.max_new_tokens = kwargs.get("max_new_tokens", 256)
        self.min_new_tokens = kwargs.get("min_new_tokens", 0)

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
                                    num_return_sequences=self.num_return_sequences, max_new_tokens=max_new_tokens, min_new_tokens=self.min_new_tokens, 
                                    pad_token_id=self.tokenizer.eos_token_id, use_cache=False)
        try:
            outputs = outputs[0][len(inputs[0]):] if not return_prompt else outputs[0]
        except:
            outputs = outputs[0][len(inputs['input_ids'][0]):] if not return_prompt else outputs[0]
        decoded_output = self.tokenizer.decode(outputs, skip_special_tokens=True)
        
        # Remove llama special words
        decoded_output = decoded_output.replace("assistant", "").replace("user", "").replace("system", "")

        return decoded_output