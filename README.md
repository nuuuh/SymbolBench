# SymbolBench

A comprehensive benchmark for evaluating Large Language Models (LLMs) on symbolic reasoning tasks across three domains: **Boolean Networks (BN)**, **Differential Equations (DE)**, and **Structural Causal Models (SCM)**.

## Overview

SymbolBench evaluates the ability of LLMs to discover and reason about symbolic expressions from data. The benchmark includes:

- **Boolean Networks (BN)**: Inference of Boolean logic rules from state transition data
- **Differential Equations (DE)**: Discovery of ODEs from time-series data (Strogatz and Physiome datasets)
- **Structural Causal Models (SCM)**: Causal discovery from observational time-series data

## Installation

```bash
# Clone the repository
git clone https://github.com/nuuuh/SymbolBench.git
cd SymbolBench

# Install dependencies
pip install -r requirements.txt
```

## Project Structure

```
SymbolBench/
├── pipeline.py              # Main pipeline for LLM-based symbolic regression
├── hybrid_pipeline.py       # Hybrid pipeline combining LLM with genetic programming
├── conf/                    # Hydra configuration files
│   ├── model/              # Model configurations (GPT, Qwen, DeepSeek, etc.)
│   ├── experiment/         # Experiment configurations (BN, DE, SCM)
│   └── dataset/            # Dataset configurations
├── scripts/                 # Experiment scripts
│   ├── BN/                 # Boolean Network experiments
│   ├── DE/                 # Differential Equation experiments
│   │   ├── strogatz/      # Strogatz dataset scripts
│   │   └── physiome/      # Physiome dataset scripts
│   └── SCM/                # Structural Causal Model experiments
├── data/                    # Datasets
├── datasets/                # Dataset loaders
├── models/                  # Model implementations
├── prompts/                 # Prompt templates for each task
├── verifiers/               # Verification modules for each task
├── analysis/                # Analysis scripts
└── summary/                 # Result summarization scripts
```

## Configuration

**Step1:** Go to the conf/config.yaml and set the "root:" to your own directory.

**Step2:** If using models from OpenAI of TogetherAI, set your API keys within scripts under the ***script*** folder:

```bash
export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"
```

## Running Experiments

### Boolean Networks (BN)

```bash
# Base setting (data only)
bash scripts/BN/BN_base.sh <start_idx> <end_idx> <batch_size>

# Context setting (data + domain labels)
bash scripts/BN/BN_context.sh <start_idx> <end_idx> <batch_size> <model>

# Reasoning setting (with chain-of-thought)
bash scripts/BN/BN_reasoning.sh <start_idx> <end_idx> <batch_size> <model>

# Hybrid approach (LLM + Genetic Programming)
bash scripts/BN/hybrid_LLM.sh <start_idx> <end_idx> <batch_size> <model>
```

### Differential Equations (DE)

**Strogatz Dataset:**
```bash
# Base setting
bash scripts/DE/strogatz/DE_base.sh <start_idx> <end_idx> <batch_size> <model>

# Context setting
bash scripts/DE/strogatz/DE_context.sh <start_idx> <end_idx> <batch_size> <model>

# Reasoning setting
bash scripts/DE/strogatz/DE_reasoning.sh <start_idx> <end_idx> <batch_size> <model>
```


### Structural Causal Models (SCM)

```bash
# Base setting
bash scripts/SCM/SCM_base.sh <start_idx> <end_idx> <batch_size>

# Context setting
bash scripts/SCM/SCM_context.sh <start_idx> <end_idx> <batch_size> <model>

# Reasoning setting
bash scripts/SCM/SCM_reasoning.sh <start_idx> <end_idx> <batch_size> <model>
```

## Supported Models

The benchmark supports various LLMs through configuration files in `conf/model/`. Currently, models from **OpenAI**, **TogetherAI**, and **HuggingFace** are supported


## Evaluation

### Evaluate Boolean Networks
```bash
python evaluate_bn.py --OOD 0
```

### Evaluate Differential Equations
```bash
# Strogatz dataset
python evaluate_DEs.py --dataset strogatz_extended --model <model_name>
```

### Evaluate Structural Causal Models
```bash
python evaluate_scm.py
```

## Experiment Settings

Each task supports four experimental settings:

| Setting | Description | Use Context | Reasoning |
|---------|-------------|-------------|-----------|
| **Naive** | Direct prediction without data | ❌ | ❌ |
| **Base** | With historical experience | ✅ | ✅ / ❌ |
| **Context** | Data + domain knowledge | ✅ | ✅ / ❌ |
| **Reasoning** | Data + domain knowledge + CoT | ✅ | ✅ |

## Datasets

- **Boolean Networks**: 65 biological regulatory networks from BoolNet
- **Differential Equations**: 
  - Strogatz: 63 ODEs from nonlinear dynamics
  - Physiome: Curated biomedical ODEs
- **Structural Causal Models**: 190 causal graphs from differential equations's variable dependencies.

## Output

Results are saved in the `runs/` directory with the following structure:
```
runs/<model>/<modality>/<experiment_name>/Expr_<idx>/
├── final_results.npy 
└── ...
```

## Contributions

We adopt part of the code structure from *["In-Context Symbolic Regression: Leveraging Large Language Models for Function Discovery"](https://arxiv.org/abs/2404.19094)*

## Citation

If you use SymbolBench in your research, please cite:

```bibtex
@article{liu2025can,
  title={Can Large Language Models Adequately Perform Symbolic Reasoning Over Time Series?},
  author={Liu, Zewen and Ni, Juntong and Tang, Xianfeng and Lau, Max SY and Yin, Wenpeng and Jin, Wei},
  journal={arXiv preprint arXiv:2508.03963},
  year={2025}
}
```

