import os
import inspect
import importlib.util
import re
import math
import json

def load_module_from_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

results = []
root = os.path.dirname(__file__)
models_root = os.path.join(root, 'models')

for domain in sorted(os.listdir(models_root)):
    domain_path = os.path.join(models_root, domain)
    if not os.path.isdir(domain_path):
        continue
    for model in sorted(os.listdir(domain_path)):
        model_path = os.path.join(domain_path, model)
        model_file = os.path.join(model_path, 'model.py')
        # initialize entry and ensure initial_states always present
        entry = {'domain': domain, 'model': model, 'initial_states': []}
        if not os.path.isfile(model_file):
            entry['error'] = 'model.py not found'
            results.append(entry)
            continue
        try:
            # load module
            module_name = f"models.{domain}.{model}.model"
            mod = load_module_from_path(module_name, model_file)
            # extract raw sources
            raw_ode = inspect.getsource(mod.computeRates)
            raw_alg = inspect.getsource(mod.computeAlgebraic)
            # parse algebraic assignments to inline in ODEs
            alg_map = {}
            for line in raw_alg.splitlines():
                m_alg = re.match(r"\s*algebraic\[(\d+)\]\s*=\s*(.+)", line)
                if m_alg:
                    alg_map[int(m_alg.group(1))] = m_alg.group(2).strip()
            # fully expand algebraic definitions recursively
            expanded_alg = {}
            for idx in sorted(alg_map):
                expr = alg_map[idx]
                def _rep_nested(m): return f"({expanded_alg[int(m.group(1))]})"
                while re.search(r"algebraic\[(\d+)\]", expr):
                    expr = re.sub(r"algebraic\[(\d+)\]", _rep_nested, expr)
                expanded_alg[idx] = expr
            # inline algebraic formulas and drop algebraic lines in ODE source
            ode_lines = []
            for line in raw_ode.splitlines():
                # skip algebraic assignments
                if re.match(r"\s*algebraic\[", line):
                    continue
                ode_lines.append(line)
            ode_text = "\n".join(ode_lines)
            # replace algebraic[...] with its expanded expression
            def _rep_alg(m):
                idx = int(m.group(1))
                return f"({expanded_alg.get(idx, m.group(0))})"
            entry['ode_source'] = re.sub(r"algebraic\[(\d+)\]", _rep_alg, ode_text)
            # transform power() to **
            entry['ode_source'] = re.sub(r"power\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1**\2)", entry['ode_source'])
            # replace relational calls
            entry['ode_source'] = re.sub(r"greater_equal\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1>=\2)", entry['ode_source'])
            entry['ode_source'] = re.sub(r"less\(\s*([^,]+?)\s*,\s*([^\)]+?)\s*\)", r"(\1<\2)", entry['ode_source'])
            # replace piecewise calls with Python ternary expressions
            entry['ode_source'] = re.sub(
                r"custom_piecewise\(\s*\[\s*(.*?)\s*,\s*(.*?)\s*,\s*True\s*,\s*(.*?)\s*\]\)",
                r"(\2 if \1 else \3)", entry['ode_source'])
            # replace bitwise AND with and
            entry['ode_source'] = entry['ode_source'].replace('&', ' and ')
            # drop separate algebraic_source since inlined
            entry.pop('algebraic_source', None)
            # init constants and states with fallback for ArrayType errors
            init_src = inspect.getsource(mod.initConsts)
            try:
                flag=True
                # if 'yang_tong' in model:
                #     import ipdb; ipdb.set_trace()
                # if 'ostby_oyehaug_einevoll' in model:
                #         import ipdb; ipdb.set_trace()

                states, constants = mod.initConsts()
                # Convert numpy arrays and numpy scalars to native Python lists of floats
                raw_states = states.tolist() if hasattr(states, 'tolist') else states
                initial_states = [float(v) for v in raw_states]
                raw_constants = constants.tolist() if hasattr(constants, 'tolist') else constants
                constants_list = [float(v) for v in raw_constants]
                # detect swapped tuple: if initial_states length matches constants count
                sizeStates_attr = getattr(mod, 'sizeStates', None)
                sizeStates = sizeStates_attr if sizeStates_attr is not None else len(entry.get('legend_states', []))
                if len(initial_states) != sizeStates and len(constants_list) == sizeStates:
                    # swap lists
                    entry['initial_states'] = constants_list
                    entry['constants'] = initial_states
                else:
                    entry['initial_states'] = initial_states
                    entry['constants'] = constants_list
            except TypeError as e:
                flag=False

                sizeStates = getattr(mod, 'sizeStates', len(entry.get('legend_states', [])))
                sizeConsts = getattr(mod, 'sizeConstants', len(entry.get('legend_constants', [])))
                initial_states = [0.0] * sizeStates
                constants_list = [0.0] * sizeConsts
                
                for m in re.finditer(r'states\[(\d+)\]\s*=\s*([0-9.eE+-]+)', init_src):
                    idx, val = int(m.group(1)), float(m.group(2))
                    if idx < sizeStates:
                        initial_states[idx] = val
                # parse constant assignments (including expressions)
                for m in re.finditer(r'constants\[(\d+)\]\s*=\s*([^\n]+)', init_src):
                    idx = int(m.group(1))
                    expr = m.group(2).strip()
                    if idx < sizeConsts:
                        try:
                            val = eval(expr, {'__builtins__':None, 'log':math.log, 'pow':math.pow}, {'constants': constants_list})
                            constants_list[idx] = float(val)
                        except Exception:
                            pass
                # import ipdb; ipdb.set_trace()
                entry['initial_states'] = initial_states
                entry['constants'] = constants_list
                flag=True

            if flag == False:
                print(f"exception not resolved for {model}")
            # ensure initial_states list is not empty
            assert 'initial_states' in entry, "initial_states not found"
            # validate initial_states length against module sizeStates or legend_states
            sizeStates = getattr(mod, 'sizeStates', None)
            if sizeStates is None:
                sizeStates = len(entry.get('legend_states', []))
            if len(entry['initial_states']) != sizeStates:
                entry['init_error'] = (
                    f"initial_states length mismatch: expected {sizeStates}, got {len(entry['initial_states'])}"
                )
            # extract timespan from solve_model source
            solve_src = inspect.getsource(mod.solve_model)
            entry['solve_model_source'] = solve_src
            # find timespan expression
            for line in solve_src.splitlines():
                if 'voi' in line and '=' in line:
                    part = line.strip().split('=', 1)[1].strip()
                    if part.startswith('linspace') or 'array' in part:
                        entry['timespan'] = part
                        break
            # extract legends describing states, algebraic variables, voi and constants
            try:
                ls, la, lv, lc = mod.createLegends()
                entry['legend_states'] = ls
                entry['legend_algebraic'] = la
                entry['legend_voi'] = lv
                entry['legend_constants'] = lc
            except Exception as e:
                entry['legend_error'] = str(e)
        except Exception as e:
            entry['error'] = str(e)
        if 'error' in entry:
            print(f"[DEBUG] Error processing {domain}.{model}: {entry['error']}")
            continue
        results.append(entry)

# save output
out_file = os.path.join(root, 'extracted_odes.json')
with open(out_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Extraction complete: {len(results)} models processed. Output written to {out_file}")
