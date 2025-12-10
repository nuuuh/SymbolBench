#!/usr/bin/env python
import json
import re


def transform_source(source: str, n_constants: int) -> str:
    """
    Extract only rates assignments, order by variable index, and replace states[], constants[], and algebraic[] with generic symbols.
    n_constants: number of constant symbols before algebraic terms
    """
    lines = source.splitlines()
    # extract computed constant assignments
    const_lines = []
    for line in lines:
        m_const = re.match(r"\s*constants\[(\d+)\]\s*=\s*(.+)", line)
        if m_const:
            idx_c = int(m_const.group(1))
            expr_c = m_const.group(2).strip()
            const_lines.append((idx_c, expr_c))
    # sort computed constants by index
    const_lines.sort(key=lambda x: x[0])
    # transform computed constant expressions
    transformed_consts = []
    for idx_c, expr_c in const_lines:
        expr_c = re.sub(r"states\[(\d+)\]", lambda m: f"x_{int(m.group(1))}", expr_c)
        expr_c = re.sub(r"constants\[(\d+)\]", lambda m: f"c_{int(m.group(1))}", expr_c)
        expr_c = re.sub(r"algebraic\[(\d+)\]", lambda m: f"c_{n_constants+int(m.group(1))}", expr_c)
        expr_c = re.sub(r"power\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1**\2)", expr_c)
        # replace relational functions with operators
        expr_c = re.sub(r"greater_equal\(\s*([^,]+?)\s*,\s*([^\)]+?)\)", r"(\1>=\2)", expr_c)
        expr_c = re.sub(r"less\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1<\2)", expr_c)
        # replace piecewise call with Python conditional
        expr_c = re.sub(
            r"custom_piecewise\(\s*\[\s*([^,]+?)\s*,\s*([^,]+?)\s*,\s*True\s*,\s*([^\]]+?)\s*\]\)",
            r"(\2 if \1 else \3)", expr_c)
        # replace bitwise AND in conditions with Python and
        expr_c = expr_c.replace('&', ' and ')
        transformed_consts.append(expr_c)

    # extract computed algebraic assignments
    alg_lines = []
    for line in lines:
        m_alg = re.match(r"\s*algebraic\[(\d+)\]\s*=\s*(.+)", line)
        if m_alg:
            idx_a = int(m_alg.group(1))
            expr_a = m_alg.group(2).strip()
            alg_lines.append((idx_a, expr_a))
    alg_lines.sort(key=lambda x: x[0])
    transformed_algs = []
    alg_exprs = {}
    # if model_name == 'grange_2001':
    #     import ipdb; ipdb.set_trace()
    for idx_a, expr_a in alg_lines:
        expr_a_copy = expr_a
        # replace states and constants
        expr = re.sub(r"states\[(\d+)\]", lambda m: f"x_{int(m.group(1))}", expr_a_copy)
        expr = re.sub(r"constants\[(\d+)\]", lambda m: f"c_{int(m.group(1))}", expr)
        # recursively expand algebraic references
        def _rep(m): return f"({alg_exprs[int(m.group(1))]})"
        while re.search(r"algebraic\[(\d+)\]", expr):
            expr = re.sub(r"algebraic\[(\d+)\]", _rep, expr)
        # apply other transformations
        expr = re.sub(r"power\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1**\2)", expr)
        expr = re.sub(r"greater_equal\(\s*([^,]+?)\s*,\s*([^\)]+?)\)", r"(\1>=\2)", expr)
        expr = re.sub(r"less\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1<\2)", expr)
        expr = re.sub(
            r"custom_piecewise\(\s*\[\s*([^,]+?)\s*,\s*([^,]+?)\s*,\s*True\s*,\s*([^\]]+?)\s*\]\)",
            r"(\2 if \1 else \3)", expr)
        expr = expr.replace('&', ' and ')
        alg_exprs[idx_a] = expr
        transformed_algs.append(expr)

    # collect rate expressions
    rate_lines = []
    for line in lines:
        m = re.match(r"\s*rates\[(\d+)\]\s*=\s*(.+)", line)
        if m:
            idx = int(m.group(1))
            expr = m.group(2).strip()
            rate_lines.append((idx, expr))
    # sort by state index
    rate_lines.sort(key=lambda x: x[0])
    # apply replacements and build ordered equations
    transformed = []
    # if model_name == 'grange_2001':
    #     import ipdb; ipdb.set_trace()
    for idx, expr in rate_lines:
        expr = re.sub(r"states\[(\d+)\]", lambda m: f"x_{int(m.group(1))}", expr)
        expr = re.sub(r"constants\[(\d+)\]", lambda m: f"c_{int(m.group(1))}", expr)
        # remove direct algebraic mapping since we expanded them
        # other transformations...
        expr = re.sub(r"power\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1**\2)", expr)
        expr = re.sub(r"greater_equal\(\s*([^,]+?)\s*,\s*([^\)]+?)\)", r"(\1>=\2)", expr)
        expr = re.sub(r"less\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1<\2)", expr)
        expr = re.sub(
            r"custom_piecewise\(\s*\[\s*([^,]+?)\s*,\s*([^,]+?)\s*,\s*True\s*,\s*([^\]]+?)\s*\]\)",
            r"(\2 if \1 else \3)", expr)
        expr = expr.replace('&', ' and ')
        transformed.append(expr)
    # combine computed constant, inlined algebraic, and rate equations
    combined = "\n".join(transformed_consts + transformed_algs + transformed)
    # final ensure all power() calls are converted to **
    combined = re.sub(r"power\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1**\2)", combined)
    return combined


def main():
    with open('extracted_odes.json') as f:
        data = json.load(f)


    idx = 0

    generic = []
    for entry in data:
        if 'ode_source' not in entry:
            continue
        name = f"{entry['domain']}.{entry['model']}"
        ode_src = entry['ode_source']

        global model_name

        model_name = entry['model']

        try:
            n_consts = len(entry.get('legend_constants', []))
            transformed = transform_source(ode_src, n_consts)
        except Exception as e:
            print(f"[DEBUG] Error transforming {name}: {e}")
            continue
        # build mappings from new symbols to legend descriptions
        # import ipdb; ipdb.set_trace()
        x_map = {f"x_{i}": desc for i, desc in enumerate(entry.get('legend_states', []) )}
        # map constants and algebraic terms into unified c_map
        c_map = {}
        consts = entry.get('legend_constants', [])
        algs = entry.get('legend_algebraic', [])
        for j, desc in enumerate(consts):
            c_map[f"c_{j}"] = desc
        for k, desc in enumerate(algs):
            c_map[f"c_{len(consts)+k}"] = desc
        # include all original fields except ode and algebraic sources and solve_model
        entry_copy = {k: v for k, v in entry.items() if k not in ('ode_source','algebraic_source','solve_model_source')}
        # split transformed string into individual ODEs and join with '|'
        eqs = transformed.split("\n") if transformed else []
        # assert number of equations matches number of state variables
        expected = len(entry.get('initial_states', []))

        assert len(x_map) == len(eqs), (
                f"Dim mismatch for {name}: expected {expected} equations, got {len(eqs)}"
            )

        assert expected == len(eqs), (
            f"Dim mismatch for {name}: expected {expected} equations, got {len(eqs)}"
        )
        entry_copy['id'] = idx
        idx += 1
        entry_copy['equation'] = "|".join(eqs)
        # record number of equations
        entry_copy['dim'] = len(eqs)
        entry_copy['x_mapping'] = x_map
        entry_copy['c_mapping'] = c_map
        generic.append(entry_copy)

    # write out generic ODE definitions
    with open('generic_odes.json', 'w') as f:
        json.dump(generic, f, indent=2)

    print(f"Transformed {len(generic)} ODEs into generic_odes.json")

    print("num dim=1: ", len([entry for entry in generic if entry.get('dim', 0) == 1]))
    print("num dim=2: ",len([entry for entry in generic if entry.get('dim', 0) == 2]))
    print("num dim=3: ",len([entry for entry in generic if entry.get('dim', 0) == 3]))
    print("num dim=4: ",len([entry for entry in generic if entry.get('dim', 0) == 4]))

    # extract and save samples with dim <= 4 and no logical operators
    small = [entry for entry in generic
             if entry.get('dim', 0) <= 4
             and entry.get('dim', 0) > 0
             and '>' not in entry.get('equation', '')
             and '<' not in entry.get('equation', '')
             and 'voi' not in entry.get('equation', '')
             and 'and' not in entry.get('equation', '')]
    small = sorted(small, key=lambda x: x['dim'])

    with open('small_odes.json', 'w') as f_small:
        json.dump(small, f_small, indent=2)
    print(f"Saved {len(small)} small ODEs into small_odes.json")


if __name__ == '__main__':
    main()
