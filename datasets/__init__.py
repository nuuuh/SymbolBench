from .DE_dataset import DE_dataset
from .BN_dataset import BN_dataset
from .BN_dataset_0 import BN_dataset_0
from .SCM_dataset import SCM_dataset
from .data_wrapper import DataWrapper


def get_dataset(cfg, logger):
    if cfg.dataset.name == "DE_dataset":
        dataset = DE_dataset(cfg, logger)
    elif cfg.dataset.name == "BN_dataset":
        dataset = BN_dataset(cfg, logger)
    elif cfg.dataset.name == "SCM_dataset":
        dataset = SCM_dataset(cfg, logger)
    else:
        raise ValueError(f"Unknown dataset: {cfg['dataset']}")
    
    return dataset