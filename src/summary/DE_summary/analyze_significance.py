#!/usr/bin/env python3
"""
Analyze statistical significance results from DE t-test files.
Only mark values with p < 0.05 as significant.
"""

import csv
import os
import glob
from collections import defaultdict

def parse_significance_files():
    """Parse all significant test files and extract p < 0.05 results."""
    
    # Find all significant test files
    pattern = "summary/DE_summary/DE_ttest_*_significant_tests.csv"
    files = glob.glob(pattern)
    
    print(f"Found {len(files)} significance test files:")
    for f in files:
        print(f"  {f}")
    
    # Dictionary to store results by model
    results = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))
    
    for file_path in files:
        # Extract model name from filename
        filename = os.path.basename(file_path)
        # Format: DE_ttest_{model}_all_vs_pysr_textual_significant_tests.csv
        model_name = filename.replace('DE_ttest_', '').replace('_all_vs_pysr_textual_significant_tests.csv', '')
        
        print(f"\nProcessing {model_name}...")
        
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            
            print(f"  Found {len(rows)} total significant tests")
            
            # Filter for p < 0.05
            significant_rows = [row for row in rows if float(row['p_value']) < 0.05]
            print(f"  Found {len(significant_rows)} tests with p < 0.05")
            
            for row in significant_rows:
                split = row['split']
                metric = row['metric']
                dimension = int(row['dimension'])
                eval_type = row['source_eval_type']
                p_value = float(row['p_value'])
                
                # Store the evaluation type for this significant result
                results[model_name][split][metric][dimension][eval_type] = p_value
                
        except Exception as e:
            print(f"  Error processing {file_path}: {e}")
    
    return results

def create_summary_table(results):
    """Create a summary table showing which combinations have p < 0.05."""
    
    models = sorted(results.keys())
    splits = ['ID', 'OOD']
    metrics = ['sr2', 'ACC_0.9']
    dimensions = [1, 2, 3, 4]
    
    print("\n" + "="*100)
    print("STATISTICAL SIGNIFICANCE SUMMARY (p < 0.05)")
    print("="*100)
    
    # Create markdown table
    markdown_lines = [
        "# DE T-Test Significance Summary (p < 0.05)",
        "",
        "## Statistical Significance Results (p < 0.05 only)",
        "",
        "| Model | Split | Metric | Dim 1 | Dim 2 | Dim 3 | Dim 4 | Evaluation Types with Significance |",
        "|-------|-------|--------|-------|-------|-------|-------|-----------------------------------|"
    ]
    
    for model in models:
        model_results = results[model]
        
        for split in splits:
            for metric in metrics:
                dim_results = []
                all_eval_types = set()
                
                for dim in dimensions:
                    if (split in model_results and 
                        metric in model_results[split] and 
                        dim in model_results[split][metric]):
                        
                        eval_types = list(model_results[split][metric][dim].keys())
                        if eval_types:
                            dim_results.append("✓")
                            all_eval_types.update(eval_types)
                        else:
                            dim_results.append("-")
                    else:
                        dim_results.append("-")
                
                # Format evaluation types
                eval_types_str = ", ".join(sorted(all_eval_types)) if all_eval_types else "-"
                
                # Only add row if there's at least one significant result
                if any(r == "✓" for r in dim_results):
                    if split == splits[0] and metric == metrics[0]:  # First row for model
                        model_display = f"**{model}**"
                    else:
                        model_display = ""
                    
                    row = f"| {model_display} | {split} | {metric} | {' | '.join(dim_results)} | {eval_types_str} |"
                    markdown_lines.append(row)
    
    # Add detailed statistics
    markdown_lines.extend([
        "",
        "## Detailed Statistics:",
        ""
    ])
    
    for model in models:
        model_results = results[model]
        markdown_lines.append(f"### {model}:")
        
        total_significant = 0
        for split in splits:
            for metric in metrics:
                for dim in dimensions:
                    if (split in model_results and 
                        metric in model_results[split] and 
                        dim in model_results[split][metric]):
                        total_significant += len(model_results[split][metric][dim])
        
        markdown_lines.append(f"- Total significant tests (p < 0.05): {total_significant}")
        
        # Show breakdown by split/metric
        for split in splits:
            for metric in metrics:
                count = 0
                dims_with_sig = []
                for dim in dimensions:
                    if (split in model_results and 
                        metric in model_results[split] and 
                        dim in model_results[split][metric]):
                        count += len(model_results[split][metric][dim])
                        dims_with_sig.append(str(dim))
                
                if count > 0:
                    markdown_lines.append(f"  - {split} {metric}: {count} significant tests (dims: {', '.join(dims_with_sig)})")
        
        markdown_lines.append("")
    
    # Add summary patterns
    markdown_lines.extend([
        "## Summary Patterns (p < 0.05):",
        "",
        "### By Model:",
    ])
    
    for model in models:
        model_results = results[model]
        total_sig = 0
        for split in model_results:
            for metric in model_results[split]:
                for dim in model_results[split][metric]:
                    total_sig += len(model_results[split][metric][dim])
        markdown_lines.append(f"- **{model}**: {total_sig} significant tests")
    
    return "\n".join(markdown_lines)

def main():
    """Main function to analyze significance and create summary."""
    
    print("Analyzing DE t-test significance results...")
    print("Significance threshold: p < 0.05")
    
    # Parse all significance files
    results = parse_significance_files()
    
    if not results:
        print("No significant results found!")
        return
    
    # Create summary table
    summary_markdown = create_summary_table(results)
    
    # Save to file
    output_file = "summary/DE_summary/t_test_significance_summary_p005.md"
    with open(output_file, 'w') as f:
        f.write(summary_markdown)
    
    print(f"\nSummary saved to: {output_file}")
    print("\nPreview of results:")
    print("="*50)
    print(summary_markdown[:2000] + "..." if len(summary_markdown) > 2000 else summary_markdown)

if __name__ == "__main__":
    main()
