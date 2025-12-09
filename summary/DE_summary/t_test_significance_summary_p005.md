# DE T-Test Significance Summary (p < 0.05)

## Statistical Significance Results (p < 0.05 only)

| Model | Split | Metric | Dim 1 | Dim 2 | Dim 3 | Dim 4 | Evaluation Types with Significance |
|-------|-------|--------|-------|-------|-------|-------|-----------------------------------|
| **ChatTS-14B** | ID | sr2 | ✓ | ✓ | ✓ | ✓ | context, reasoning |
|  | OOD | sr2 | ✓ | - | - | - | context, reasoning |
|  | OOD | ACC_0.9 | ✓ | - | - | - | context, reasoning |
| **Llama-3.2-3B** | ID | sr2 | ✓ | ✓ | - | ✓ | context, reasoning |
|  | ID | ACC_0.9 | - | - | ✓ | ✓ | context, reasoning |
|  | OOD | sr2 | ✓ | - | ✓ | - | context, reasoning |
|  | OOD | ACC_0.9 | ✓ | - | - | - | context, reasoning |
| **Mathstral-7B-v0.1** | ID | sr2 | ✓ | ✓ | ✓ | - | context, reasoning |
|  | OOD | sr2 | ✓ | - | - | - | context, reasoning |
|  | OOD | ACC_0.9 | ✓ | - | - | - | context, reasoning |
| **Qwen2.5-14B** | ID | sr2 | ✓ | ✓ | ✓ | ✓ | base, context, naive, reasoning |
|  | ID | ACC_0.9 | - | - | - | ✓ | base, context, naive, reasoning |
|  | OOD | sr2 | ✓ | ✓ | ✓ | - | base, context, naive, reasoning |
|  | OOD | ACC_0.9 | ✓ | ✓ | ✓ | - | base, context, naive, reasoning |
| **gpt-4o-mini** | ID | sr2 | ✓ | ✓ | ✓ | ✓ | context, reasoning |
|  | ID | ACC_0.9 | - | ✓ | - | - | context, reasoning |
|  | OOD | sr2 | ✓ | ✓ | - | - | context, reasoning |
|  | OOD | ACC_0.9 | ✓ | ✓ | - | - | context, reasoning |

## Detailed Statistics:

### ChatTS-14B:
- Total significant tests (p < 0.05): 10
  - ID sr2: 6 significant tests (dims: 1, 2, 3, 4)
  - OOD sr2: 2 significant tests (dims: 1)
  - OOD ACC_0.9: 2 significant tests (dims: 1)

### Llama-3.2-3B:
- Total significant tests (p < 0.05): 12
  - ID sr2: 5 significant tests (dims: 1, 2, 4)
  - ID ACC_0.9: 2 significant tests (dims: 3, 4)
  - OOD sr2: 3 significant tests (dims: 1, 3)
  - OOD ACC_0.9: 2 significant tests (dims: 1)

### Mathstral-7B-v0.1:
- Total significant tests (p < 0.05): 9
  - ID sr2: 5 significant tests (dims: 1, 2, 3)
  - OOD sr2: 2 significant tests (dims: 1)
  - OOD ACC_0.9: 2 significant tests (dims: 1)

### Qwen2.5-14B:
- Total significant tests (p < 0.05): 29
  - ID sr2: 13 significant tests (dims: 1, 2, 3, 4)
  - ID ACC_0.9: 4 significant tests (dims: 4)
  - OOD sr2: 6 significant tests (dims: 1, 2, 3)
  - OOD ACC_0.9: 6 significant tests (dims: 1, 2, 3)

### gpt-4o-mini:
- Total significant tests (p < 0.05): 17
  - ID sr2: 8 significant tests (dims: 1, 2, 3, 4)
  - ID ACC_0.9: 2 significant tests (dims: 2)
  - OOD sr2: 4 significant tests (dims: 1, 2)
  - OOD ACC_0.9: 3 significant tests (dims: 1, 2)

## Summary Patterns (p < 0.05):

### By Model:
- **ChatTS-14B**: 10 significant tests
- **Llama-3.2-3B**: 12 significant tests
- **Mathstral-7B-v0.1**: 9 significant tests
- **Qwen2.5-14B**: 29 significant tests
- **gpt-4o-mini**: 17 significant tests