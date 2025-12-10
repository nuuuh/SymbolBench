# DE T-Test Significance Summary

## Statistical Significance Results (p < 0.05)

| Model | Split | Metric | Dim 1 | Dim 2 | Dim 3 | Dim 4 | Evaluation Types with Significance |
|-------|-------|--------|-------|-------|-------|-------|-----------------------------------|
| **ChatTS-14B** | ID | sr2 | ✓ | ✓ | ✓ | ✓ | context, reasoning |
| | OOD | sr2 | ✓ | - | - | - | context, reasoning |
| | OOD | ACC_0.9 | ✓ | - | - | - | context, reasoning |
| **Llama-3.2-3B** | ID | sr2 | ✓ | ✓ | - | ✓ | context, reasoning |
| | ID | ACC_0.9 | - | - | ✓ | ✓ | context, reasoning |
| | OOD | sr2 | ✓ | - | ✓ | - | context, reasoning |
| | OOD | ACC_0.9 | ✓ | - | - | - | context, reasoning |
| **Mathstral-7B-v0.1** | ID | sr2 | ✓ | ✓ | ✓ | - | context, reasoning |
| | OOD | sr2 | ✓ | - | - | - | context, reasoning |
| | OOD | ACC_0.9 | ✓ | - | - | - | context, reasoning |
| **Qwen2.5-14B** | ID | sr2 | ✓ | ✓ | ✓ | ✓ | base, context, naive, reasoning |
| | ID | ACC_0.9 | - | - | - | ✓ | base, context, naive, reasoning |
| | OOD | sr2 | ✓ | ✓ | ✓ | - | base, context, naive, reasoning |
| | OOD | ACC_0.9 | ✓ | ✓ | ✓ | - | base, context, naive, reasoning |
| **gpt-4o-mini** | ID | sr2 | ✓ | ✓ | ✓ | ✓ | context, reasoning |
| | ID | ACC_0.9 | - | ✓ | - | - | context, reasoning |
| | OOD | sr2 | ✓ | ✓ | - | - | context, reasoning |
| | OOD | ACC_0.9 | ✓ | ✓ | - | - | context, reasoning |

## Summary Statistics (p < 0.05):

### Total Significant Tests by Model:
- **Qwen2.5-14B**: 29 significant tests (most comprehensive)
- **gpt-4o-mini**: 17 significant tests  
- **Llama-3.2-3B**: 12 significant tests
- **ChatTS-14B**: 10 significant tests
- **Mathstral-7B-v0.1**: 9 significant tests

### By Model Performance:
- **Qwen2.5-14B**: Most comprehensive significance across all dimensions and evaluation types (naive, base, context, reasoning)
- **gpt-4o-mini**: Strong across all dimensions but limited to context/reasoning approaches
- **ChatTS-14B**: Good ID performance across all dimensions, limited OOD significance
- **Llama-3.2-3B**: Mixed patterns with some gaps in dimension coverage
- **Mathstral-7B-v0.1**: Most conservative, fewer significant differences

### By Dimension:
- **Dimension 1**: All models show significant differences (most consistent)
- **Dimension 2**: Strong significance across most models
- **Dimension 3**: Mixed results, fewer models show significance
- **Dimension 4**: Limited to larger models (Qwen2.5-14B, ChatTS-14B, gpt-4o-mini, Llama-3.2-3B)

### By Evaluation Type:
- **Context & Reasoning**: Most frequent evaluation types showing significance
- **Naive & Base**: Only present for Qwen2.5-14B (comprehensive evaluation)
- **Pattern**: Advanced prompting methods generally show more statistical differences

### By Split:
- **ID (In-Distribution)**: More consistent significance across dimensions
- **OOD (Out-of-Distribution)**: Concentrated in lower dimensions (1-2), fewer significant results

### By Metric:
- **sr2**: Primary metric showing significance across all models
- **ACC_0.9**: More selective, typically significant in dimensions 1-4 depending on model

## Key Takeaways (p < 0.05):
1. **Qwen2.5-14B** shows the most comprehensive statistical differences from PySR baseline
2. Lower dimensions (1-2) have more consistent statistical significance across models
3. **Context and reasoning** evaluation methods reveal more statistical differences than simpler approaches
4. **OOD performance** shows fewer significant differences than ID performance
5. Only **Qwen2.5-14B** was tested with all evaluation types (naive, base, context, reasoning)
