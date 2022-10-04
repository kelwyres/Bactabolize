import json
import sys
import os

import cobra.io
from cobra.io import read_sbml_model


from . import media_definitions
from . import util


def run(config):
    # pylint: disable=too-many-branches,too-many-statements
    print('\n========================================')
    print('Patching model ' + os.path.splitext(os.path.basename(config.draft_model_fp))[0])
    print('========================================')
    # Read in models and patch file
    with config.draft_model_fp.open('r') as fh:
        if config.draft_model_fp.suffix == '.json':
            model_draft = cobra.io.load_json_model(fh)
        elif config.draft_model_fp.suffix == '.xml':
            model_draft = read_sbml_model(fh)

    with config.ref_model_fp.open('r') as fh:
        if config.ref_model_fp.suffix == '.json':
            model_ref = cobra.io.load_json_model(fh)
        elif config.ref_model_fp.suffix == '.xml':
            model_ref = read_sbml_model(fh)

    patch = parse_patch(config.patch_fp, model_draft.id)
    # Apply patch
    for reaction_id, op in patch['reactions'].items():
        if op == 'remove':
            model_draft.remove_reactions([reaction_id])
        elif op == 'add':
            reaction = model_ref.reactions.get_by_id(reaction_id)
            model_draft.add_reaction(reaction)
        else:
            print(f'error: got bad operation {op} for {reaction_id}', file=sys.stderr)
            sys.exit(1)
    #    biomass_reaction = model_draft.reactions.get_by_id(config.biomass_reaction_id)
    #    for metabolite_id, op in patch['biomass_metabolites'].items():
    #        if op == 'remove':
    #            metabolite = model_draft.metabolites.get_by_id(metabolite_id)
    #            coefficient = biomass_reaction.metabolites[metabolite]
    #            biomass_reaction.subtract_metabolites({metabolite: coefficient})
    #        elif op == 'add':
    #            # NOTE: need to check if we'd ever need to add a biomass metabolite
    #            raise NotImplementedError
    #        else:
    #            print(f'error: got bad operation {op} for {metabolite_id}', file=sys.stderr)
    #            sys.exit(1)
    # Write model to disk
    with config.output_fp.open('w') as fh:
        cobra.io.save_json_model(model_draft, fh)
        cobra.io.write_sbml_model(model_draft, str(config.output_fp).rsplit('.', 1)[0] + '.xml')  # .xml output
    # Check if model now optimises on set media
    for reaction in model_draft.exchanges:
        reaction.lower_bound = 0
    media = media_definitions.get(config.media_type)
    for reaction_id, lower_bound in media['exchanges'].items():
        try:
            reaction = model_draft.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: draft model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound

    # Set atmospheric conditions
    reaction_oxygen = model_draft.reactions.get_by_id('EX_o2_e')
    if config.atmosphere_type == 'aerobic':
        reaction_oxygen.lower_bound = -20
    elif config.atmosphere_type == 'anaerobic':
        reaction_oxygen.lower_bound = 0
    elif config.atmosphere_type is not None:
        raise ValueError

    solution = model_draft.optimize()

    # Threshold for whether a model produces biomass
    if solution.objective_value < 1e-4:
        print(f'error: {model_draft} model failed to produce biomass on minimal media', file=sys.stderr)
        sys.exit(101)
    else:
        print(f'{model_draft} model produces biomass on minimal media')
        msg = (
            'Please cite:\n'
            '  - bioRxiv and later, published paper here\n'
            '  - Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. '
            'COBRApy: COnstraints-Based Reconstruction and Analysis for Python. '
            'BMC Syst Biol 7, 74 (2013). https://doi.org/10.1186/1752-0509-7-74'
        )
        print(msg)

    # Generate MEMOTE report file if requested
    if config.memote_report_fp:
        util.generate_memote_report(model_draft, config.memote_report_fp)


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
