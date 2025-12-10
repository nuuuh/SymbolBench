import os
import sys
import pandas as pd


class BaseVerifier:
    """
    Base class for all verifiers.
    """

    def __init__(self, cfg):
        pass

        # self.time_series
        # self.candiddate_expressions = pd.DataFrame()
        # self.additional_prompt=None
        # self.OOD = None

    def parse_model_output(self, model_output:str):
        pass

    def prepare_prompt(self, base_prompt:str)-> str:
        pass

    
    def update_candidates(self, candidates:list):
        pass

    def check_tolerance()->bool:
        pass

