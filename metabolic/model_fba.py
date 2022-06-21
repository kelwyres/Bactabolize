import itertools
import json
import sys


import cobra.io
from cobra.io import read_sbml_model


from . import fba


def run(model_fp, fba_open_value, spec_fp, output_fp):
    # pylint: disable=too-many-branches
    print('\n========================================')
    print('running FBA')
    print('========================================')
    # Read in model and spec
    with model_fp.open('r') as fh:
        if model_fp.suffix == '.json':
            model = cobra.io.load_json_model(fh)
        elif model_fp.suffix == '.xml':
            model = read_sbml_model(fh)
    spec = parse_spec(spec_fp)
    # Run FBA
    results = dict()
    for fba_name, fba_spec in spec.items():
        results[fba_name] = dict()
        for fba_type in fba_spec['fba_type']:
            if fba_type == 'potential_element_sources':
                fba_output = fba_potential_sources(model.copy(), fba_open_value, fba_spec)
            elif fba_type == 'defined_exchanges_only':
                fba_output = fba_media(model.copy(), fba_spec)
            else:
                assert False
            results[fba_name][fba_type] = fba_output
    # Write results
    with output_fp.open('w') as fh:
        header_tokens = ('fba_type', 'spec_name', 'atmosphere', 'exchange', 'categories', 'objective_value')
        print(*header_tokens, sep='\t', file=fh)
        for fba_name, fba_data in results.items():
            for fba_type, fba_output in fba_data.items():
                for data in fba_output:
                    if fba_type == 'potential_element_sources':
                        reaction_name, cats, atmosphere, value = data
                        print(fba_type, fba_name, atmosphere, reaction_name, cats, value, sep='\t', file=fh)
                    elif fba_type == 'defined_exchanges_only':
                        atmosphere, value = data
                        print(fba_type, fba_name, atmosphere, '-', '-', value, sep='\t', file=fh)


def fba_potential_sources(model, fba_open_value, spec):
    # pylint: disable=too-many-branches
    # For each potential source, run FBA
    fba_results = list()
    fba_data = fba.prepare_element_source_data(model)
    for reaction_name, categories in fba_data.items():
        # Run FBA for all category combinations
        for i in range(len(categories)):
            for fba_categories in itertools.combinations(categories, i + 1):
                # Close all exchanges then open media-defined exchanges
                reaction_bounds = {r.id: 0 for r in model.exchanges}
                for reaction_id, lower_bound in spec['exchanges'].items():
                    reaction_bounds[reaction_id] = lower_bound
                # Close default sources
                for source_name in fba_categories:
                    reaction_id = spec['default_element_sources'][source_name]
                    reaction_bounds[reaction_id] = 0
                # Open reaction investigated
                reaction_bounds[reaction_name] = fba_open_value
                # Run in aerobic and anaerobic atmosphere
                for atmosphere in ('aerobic', 'anaerobic'):
                    if atmosphere == 'aerobic':
                        reaction_bounds['EX_o2_e'] = -20
                    else:
                        reaction_bounds['EX_o2_e'] = 0
                    # Run FBA
                    objective_value = fba.run_fba(model, reaction_bounds)
                    data = (reaction_name, ','.join(fba_categories), atmosphere, objective_value)
                    fba_results.append(data)
    return fba_results


def fba_media(model, spec):
    fba_results = list()
    # Close all exchanges then open media-defined exchanges
    reaction_bounds = {r.id: 0 for r in model.exchanges}
    for reaction_id, lower_bound in spec['exchanges'].items():
        reaction_bounds[reaction_id] = lower_bound
    # Aerobic and anerobic
    for atmosphere in ('aerobic', 'anaerobic'):
        if atmosphere == 'aerobic':
            reaction_bounds['EX_o2_e'] = -20
        else:
            reaction_bounds['EX_o2_e'] = 0
        # Run FBA
        objective_value = fba.run_fba(model, reaction_bounds)
        data = (atmosphere, objective_value)
        fba_results.append(data)
    return fba_results


def parse_spec(spec_fp):
    with spec_fp.open('r') as fh:
        spec = json.load(fh)
    for fba_spec in spec.values():
        validate_spec(fba_spec)
    return spec


def validate_spec(fba_spec):
    # pylint: disable=too-many-branches,too-many-statements
    # Is a dict and has required data
    if not isinstance(fba_spec, dict):
        print('error: fba spec is not a dictionary', file=sys.stderr)
        sys.exit(1)
    # Fields
    if 'fba_type' not in fba_spec:
        print('error: no fba_types defined', file=sys.stderr)
        sys.exit(1)
    if 'exchanges' not in fba_spec:
        print('error: no exchanges defined', file=sys.stderr)
        sys.exit(1)
    if 'default_element_sources' not in fba_spec:
        print('error: no default_element_sources defined', file=sys.stderr)
        sys.exit(1)
    # Field types
    if not isinstance(fba_spec['fba_type'], list):
        print('error: fba_type field is not a list', file=sys.stderr)
        sys.exit(1)
    if not isinstance(fba_spec['exchanges'], dict):
        print('error: exchanges field is not a dictionary', file=sys.stderr)
        sys.exit(1)
    if not isinstance(fba_spec['default_element_sources'], dict):
        print('error: default_element_sources field is not a dictionary', file=sys.stderr)
        sys.exit(1)
    # Valid FBA types to run
    fba_type_valid = {'defined_exchanges_only', 'potential_element_sources'}
    for fba_type in fba_spec['fba_type']:
        if fba_type not in fba_type_valid:
            print('error: got bad fba_type: {fba_type}', file=sys.stderr)
            sys.exit(1)
    # Exchange lower bounds are numeric
    for exchange_value in fba_spec['exchanges'].values():
        if not isinstance(exchange_value, int) and not isinstance(exchange_value, float):
            print(f'error: found non-numeric exchange LB: {exchange_value}', file=sys.stderr)
            sys.exit(1)
    # All defaults are defined
    sources_valid = {'carbon', 'phosphate', 'nitrogen', 'sulfur'}
    sources_present = set(fba_spec['default_element_sources'])
    missing = sources_valid.difference(sources_present)
    undefined = sources_present.difference(sources_valid)
    if missing:
        missing_str = ', '.join(missing)
        plurality = 'sources' if len(missing) > 1 else 'source'
        print(f'error: default {plurality} missing: {missing_str}', file=sys.stderr)
        sys.exit(1)
    if undefined:
        undefined_str = ', '.join(undefined)
        plurality = 'sources' if len(undefined) > 1 else 'source'
        print(f'error: undefined default {plurality} found: {undefined_str}', file=sys.stderr)
        sys.exit(1)
    # Default exchanges present in defined exchanges
    source_exchanges_present = set(fba_spec['default_element_sources'].values())
    missing = source_exchanges_present.difference(fba_spec['exchanges'])
    if missing:
        missing_str = ', '.join(missing)
        plurality = 'sources' if len(missing) > 1 else 'source'
        print(f'error: default {plurality} not defined in exchanges: {missing_str}', file=sys.stderr)
        sys.exit(1)
    # O2 and CO2 should not be in media
    if 'EX_o2_e' in fba_spec['exchanges']:
        print('error: O2 present in exchange list', file=sys.stderr)
        sys.exit(1)
    if 'EX_co2_e' in fba_spec['exchanges']:
        print('error: CO2 present in exchange list', file=sys.stderr)
        sys.exit(1)
