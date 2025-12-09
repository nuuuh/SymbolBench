import os
import torch
# Set environment variable to handle signals properly and import juliacall before torch
os.environ['PYTHON_JULIACALL_HANDLE_SIGNALS'] = 'yes'
try:
    import juliacall
except ImportError:
    print("Warning: juliacall not available, Julia features will be disabled")
    juliacall = None

from .openai_model import OpenAIModel
from .hf_model import HuggingFaceModel
from .hybrid_DE import Hybrid_DE
from .hybrid_BN import Hybrid_BN



def load_model(model_name: str, device, dtype, cache_dir: str = None, model_args = None, exp_name=None):
    """
    Utility to load a model from the HuggingFace model hub.
    Mostly needed to deal with LLaVA models, that are not available on the model hub yet.

    Parameters
    ----------
    model_name -> the name of the model to load.
    device -> the device to load the model on.
    dtype -> the dtype to load the model with.
    cache_dir -> the cache directory to use for the model.

    Returns
    -------
    model -> the loaded model.
    """ 
    # Remove any unsupported 'proxies' from model_args
    if model_args is None:
        model_args = {}
    else:
        model_args.pop('proxies', None)
    # import ipdb; ipdb.set_trace()
    if exp_name is not None and "hybrid_GP_DE" in exp_name:
        return Hybrid_DE(model_name, device, dtype, cache_dir, **model_args)
    elif exp_name is not None and "hybrid_GP_BN" in exp_name:
        return Hybrid_BN(model_name, device, dtype, cache_dir, **model_args)
    
    elif 'gpt' in model_name or device == torch.device('cpu'):
        model = OpenAIModel(model_name, device, dtype, cache_dir, **model_args)
    else:
        model = HuggingFaceModel(model_name, device, dtype, cache_dir, **model_args)

    return model