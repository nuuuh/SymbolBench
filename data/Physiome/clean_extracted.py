import json
import re
from tqdm import tqdm


def clean_extracted(input_path, output_path=None):
    # Load generic ODE entries
    with open(input_path) as f:
        data = json.load(f)

    cleaned = []
    for t, entry in tqdm(enumerate(data)):
        expr = entry.get('equation', '')
        # fix accidental tuple syntax '(expr, 1.0/2)' to exponentiation and any leftover power() calls
        expr = re.sub(r"\(\s*(.*?)\s*,\s*1\.0/2\s*\)", r"(\1**0.5)", expr)
        expr = re.sub(r"power\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1**\2)", expr)
        if not expr  or any(op in expr for op in ('>', '<', 'and', 'voi', 'greater', 'if', 'algebra')):
            continue

        # find all constant indices used in expression
        used = sorted({int(i) for i in re.findall(r'c_(\d+)', expr)})
        # build new constant list in order
        orig_consts = entry.get('constants', [])
        try:
            new_consts = [orig_consts[i] for i in used]
        except:
            import ipdb; ipdb.set_trace()
        # mapping old index to new index
        idx_map = {old: new for new, old in enumerate(used)}
        # reindex constants in expression
        def repl(m):
            old = int(m.group(1))
            return f"c_{idx_map[old]}"
        new_expr = re.sub(r'c_(\d+)', repl, expr)
        for i in range(len(new_consts)):
            # replace c_i with c_0, c_1, ... in expression
            # import ipdb; ipdb.set_trace()
            sub_expr = re.sub(r'c_(\d+)', str(new_consts[i]), new_expr)
        # import ipdb; ipdb.set_trace()

        # assemble cleaned entry
        cleaned.append({
            'id': entry.get('id'),
            'eq': new_expr,
            'dim': entry.get('dim'),
            'consts': [new_consts],
            'init': [entry.get('initial_states')],
            'substituted': [str(sub_expr).split("|")],
            'init_constraints': entry.get('init_constraints', ''),
            'const_constraints': entry.get('const_constraints', ''),
            'eq_description': entry.get('eq_description', ''),
            'const_description': entry.get('const_description', ''),
            'var_description': entry.get('var_description', ''),
            'source': entry.get('source', '')
        })

    if output_path:
        with open(output_path, 'w') as f:
            json.dump({'equations': cleaned}, f, indent=4)
    return cleaned


if __name__ == '__main__':
    cleaned = clean_extracted(
        input_path='small_odes.json',
        output_path='aligned_small_odes.json'
    )
    print(f"Cleaned {len(cleaned)} entries.")
    # print('equations = ' + json.dumps(cleaned, indent=4))
