import json


import cobra.io


from . import media_definitions
from . import fba


def run(model_fp, fba_types, spec_fp, output_fp):
    print('\n========================================')
    print('running FBA')
    print('========================================')
    # Read in model
    with model_fp.open('r') as fh:
        model = cobra.io.load_json_model(fh)
    # Run FBA
    results = dict()
    if not fba_types or 'individual' in fba_types:
        results['individual'] = fba_individal_sources(model.copy())
    if not fba_types or 'media' in fba_types:
        results['media'] = fba_media(model.copy())
    if spec_fp and (not fba_types or 'spec' in fba_types):
        results['spec'] = fba_spec(model.copy(), spec_fp)
    # Write results
    with output_fp.open('w') as fh:
        print('fba_type', 'name/reaction', 'categories', 'objective_value', sep='\t', file=fh)
        for fba_type, data in results.items():
            if fba_type == 'individual':
                for (reaction, cats), value in data.items():
                    print(fba_type, reaction, ','.join(cats), value, sep='\t', file=fh)
            else:
                for name, value in data.items():
                    print(fba_type, name, '-', value, sep='\t', file=fh)


def fba_individal_sources(model):
    # Get formulas of metabolites involved in extracellular exchange
    metabolite_ids = list()
    for reaction in model.exchanges:
        for m in reaction.metabolites:
            [mid] = m.annotation['metanetx.chemical']
            metabolite_ids.append(mid)
    metabolite_formulas = fba.get_formulas(metabolite_ids)
    # Execute FBAs
    fba_data, exchange_sources = fba.get_individual_data(model, metabolite_formulas)
    fba_results = dict()
    for reaction_list, categories in fba_data:
        # Create dict with lower bounds for enabled exchanges
        reaction_bounds = {r: -1000 for r in reaction_list}
        objective_value = fba.run_fba(model, reaction_bounds, exchange_block=exchange_sources)
        fba_key = (reaction_list[0], tuple(categories))
        assert fba_key not in fba_results
        fba_results[fba_key] = objective_value
    return fba_results


def fba_media(model):
    # TODO: iterate all media once others have been defined
    return {'m9': fba.run_fba(model, media_definitions.m9)}


def fba_spec(model, spec_fp):
    with spec_fp.open('r') as fh:
        spec_reaction_bounds = json.load(fh)
    fba_results = dict()
    for spec, reaction_bounds in spec_reaction_bounds.items():
        objective_value = fba.run_fba(model, reaction_bounds)
        fba_results[spec] = objective_value
    return fba_results
