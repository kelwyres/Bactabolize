import json
import sys


import cobra.io


from . import media_definitions


def run(model_fp, patch_fp, output_fp):
    print('\n========================================')
    print('Patching model')
    print('========================================')
    # Read in model and patch file
    with model_fp.open('r') as fh:
        model = cobra.io.load_json_model(fh)
    patch = parse_patch(patch_fp, model.id)
    # Apply patch
    for reaction_id, op in patch['reactions'].items():
        if op == 'remove':
            model.remove_reactions([reaction_id])
        elif op == 'add':
            # TODO: write code to create Reaction objects to be added
            raise NotImplemented
        else:
            print(f'error: got bad operation {op} for {vid}', file=sys.stderr)
            sys.exit(1)
    biomass_reaction = model.reactions.get_by_id('BIOMASS_')
    for metabolite_id, op in patch['biomass_metabolites'].items():
        if op == 'remove':
            metabolite = model.metabolites.get_by_id(metabolite_id)
            coefficient = biomass_reaction.metabolites[metabolite]
            biomass_reaction.subtract_metabolites({metabolite: coefficient})
        elif op == 'add':
            # TODO: write code to create Metabolite objects to be added
            raise NotImplemented
        else:
            print(f'error: got bad operation {op} for {vid}', file=sys.stderr)
            sys.exit(1)
    # Check if model now optimises on m9
    for reaction in model.exchanges:
        reaction.lower_bound = 0
    for reaction_id, lower_bound in media_definitions.m9.items():
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: draft model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound
    solution = model.optimize()
    print(solution)


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
