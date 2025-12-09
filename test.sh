#!/bin/bash

export OPENAI_API_KEY="your_openai_api_key_here"
# export OPENAI_ORG_ID="your-org-id-here"

# python  pipeline.py experiment.symbolic_expression.name="test" experiment.textual_input="prompts/DE/DE_textual_reasoning.json" experiment.symbolic_expression.reasoning=true experiment.symbolic_expression.use_context=true experiment.symbolic_expression.idx=48 
# python  pipeline.py \
#         experiment='BN_exp' \
#         dataset='BN' \
#         experiment.symbolic_expression.name="test" \
#         experiment.symbolic_expression.data_path="data/BN.json" \
#         experiment.textual_input="prompts/BN/BN_textual_base.json" \
#         experiment.symbolic_expression.reasoning=false \
#         experiment.symbolic_expression.use_context=false \
#         experiment.symbolic_expression.idx=1 \
#         experiment.symbolic_expression.judge.enabled=false \
#         experiment.symbolic_expression.judge.gpt_model="gpt-4.1-nano" \
#         experiment.symbolic_expression.additional_prompt='none' \



dataset='data/science_sgc_total.npy' # data/science_sgc_total.npy or data/synthetic_sgc.npy
# python  pipeline.py experiment.symbolic_expression.name="test" experiment.textual_input="prompts/DE/DE_textual_reasoning.json" experiment.symbolic_expression.reasoning=true experiment.symbolic_expression.use_context=true experiment.symbolic_expression.idx=48 
python  pipeline.py \
        experiment='SCM_exp' \
        dataset='SCM' \
        experiment.symbolic_expression.name="test" \
        experiment.symbolic_expression.data_path=$dataset \
        experiment.textual_input="prompts/SCM/SCM_textual_base.json" \
        experiment.symbolic_expression.reasoning=false \
        experiment.symbolic_expression.use_context=false \
        experiment.symbolic_expression.idx=1 \
        experiment.symbolic_expression.judge.enabled=false \
        experiment.symbolic_expression.judge.gpt_model="gpt-4.1-nano" \
        experiment.symbolic_expression.additional_prompt='none' \
