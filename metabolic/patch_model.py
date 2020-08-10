import json
import sys


import cobra.io


from . import media_definitions


def run(draft_model_fp, ref_model_fp, patch_fp, output_fp):
    print('\n========================================')
    print('Patching model')
    print('========================================')
    # Read in models and patch file
    with draft_model_fp.open('r') as fh:
        draft_model = cobra.io.load_json_model(fh)
    with ref_model_fp.open('r') as fh:
        ref_model = cobra.io.load_json_model(fh)
    patch = parse_patch(patch_fp, draft_model.id)
    # Apply patch
    for reaction_id, op in patch['reactions'].items():
        if op == 'remove':
            draft_model.remove_reactions([reaction_id])
        elif op == 'add':
            reaction = ref_model.reactions.get_by_id(reaction_id)
            draft_model.add_reaction(reaction)
        else:
            print(f'error: got bad operation {op} for {vid}', file=sys.stderr)
            sys.exit(1)
    biomass_reaction = draft_model.reactions.get_by_id('BIOMASS_')
    for metabolite_id, op in patch['biomass_metabolites'].items():
        if op == 'remove':
            metabolite = draft_model.metabolites.get_by_id(metabolite_id)
            coefficient = biomass_reaction.metabolites[metabolite]
            biomass_reaction.subtract_metabolites({metabolite: coefficient})
        elif op == 'add':
            # TODO: check if we'd ever need to add a biomass metabolite
            raise NotImplemented
        else:
            print(f'error: got bad operation {op} for {vid}', file=sys.stderr)
            sys.exit(1)
    # Check if model now optimises on m9
    for reaction in draft_model.exchanges:
        reaction.lower_bound = 0
    for reaction_id, lower_bound in media_definitions.m9.items():
        try:
            reaction = draft_model.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: draft model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound
    solution = draft_model.optimize()
    # Threshold for whether a model produces biomass
    if solution.objective_value < 1e-4:
        print('error: model failed to produce biomass on minimal media', file=sys.stderr)
        sys.exit(101)
    else:
        print('model produces biomass on minimal media')


def parse_patch(patch_fp, model_name):
    # Get patch data for model
    with patch_fp.open('r') as fh:
        patch = json.load(fh)
    if model_name not in patch:
        print(f'error: could not find model {model_name} in patch data', file=sys.stderr)
        sys.exit(1)
    patch = patch[model_name]
    # Check all ops are valid
    for section in patch:
        for vid, op in patch[section].items():
            if op in {'add', 'remove'}:
                continue
            print(f'error: got bad operation {op} for {vid}', file=sys.stderr)
            sys.exit(1)
    return patch
