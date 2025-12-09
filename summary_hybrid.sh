#!/bin/bash

# GP
python summary/DE_hybrid_summary.py \
    --eval_type "hybrid_LLM_DE_wo_llm" \
    --path "hybrid_output/gpt-4o-mini/textual/"

python summary/BN_hybrid_summary.py \
    --eval_type "hybrid_LLM_BN_wo_llm" \
    --path "hybrid_output/gpt-4o-mini/textual/"


# LLM
python summary/DE_hybrid_LLM_preprocess.py \
    --eval_type "context_strogatz" \
    --path "runs/gpt-4o-mini/textual/"

python summary/DE_hybrid_summary.py \
    --eval_type "context_strogatz" \
    --path "runs/gpt-4o-mini/textual/"

python summary/BN_hybrid_summary.py \
    --eval_type "BN_context" \
    --path "runs/gpt-4o-mini/textual/"

python summary/BN_hybrid_LLM_preprocess.py \
    --eval_type "BN_context" \
    --path "runs/gpt-4o-mini/textual/"


# LLM-as-predictor+GP
python summary/DE_hybrid_summary.py \
    --eval_type "hybrid_LLM_DE_w_llm" \
    --path "hybrid_output/gpt-4o-mini/textual/"

python summary/BN_hybrid_summary.py \
    --eval_type "hybrid_LLM_BN_w_llm" \
    --path "hybrid_output/gpt-4o-mini/textual/"

# LLM-as-judge+GP with different weights
judge_weight=0.1
python summary/DE_hybrid_summary.py \
    --eval_type "hybrid_LLM_DE_w_llm_$judge_weight" \
    --path "hybrid_output/gpt-4o-mini/textual/"

python summary/BN_hybrid_summary.py \
    --eval_type "hybrid_LLM_BN_w_llm_$judge_weight" \ 
    --path "hybrid_output/gpt-4o-mini/textual/"