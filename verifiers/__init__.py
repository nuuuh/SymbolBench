from .BN_verifier import BN_Verifier
from .DE_verifier import DE_Verifier
from .SCM_verifier import SCM_Verifier


def get_verifier(cfg, dataset, logger, device=None, dtype=None):
    """
    Get the appropriate verifier class based on the name provided.

    Args:
        name (str): The name of the verifier to retrieve.

    Returns:
        BaseVerifier: An instance of the specified verifier class.
    """
    if cfg.dataset.name == "BN_dataset":
        return BN_Verifier(cfg, dataset, logger, device)
    elif cfg.dataset.name == "DE_dataset":
        return DE_Verifier(cfg, dataset, logger, device, dtype)
    elif cfg.dataset.name == "SCM_dataset":
        return SCM_Verifier(cfg, dataset, logger, device)
    else:
        raise ValueError(f"Unknown verifier name: {cfg.dataset.name}")