import sys

import cobra.flux_analysis
import cobra.io


from . import media_definitions


def run(config):
    print('\n========================================')
    print('running SGK on ' + config.model_fp.stem)
    print('========================================')
    # Load model
    with config.model_fp.open('r') as fh:
        if config.model_fp.name.endswith('.sbml'):
            model, _ = cobra.io.validate_sbml_model(fh)
        elif config.model_fp.name.endswith('.json'):
            model = cobra.io.load_json_model(fh)
        elif config.model_fp.name.endswith('.xml'):
            model = cobra.io.read_sbml_model(fh)
        else:
            assert False

    # Set growth environment
    # NOTE(SW): returning to be explicit about in-place modification
    model = set_growth_environment(model, config.media_type, config.atmosphere_type)

    # Run single gene knockout analysis and write results to disk
    results = cobra.flux_analysis.single_gene_deletion(model)
    results.to_csv(config.output_fp, sep='\t', index=False)
    print(config.model_fp.stem + ' Single Gene Knockout analysis complete')


def set_growth_environment(model, media_type, atmosphere_type):
    # Set media type
    for reaction in model.exchanges:
        reaction.lower_bound = 0
    media = media_definitions.get(media_type)
    for reaction_id, lower_bound in media['exchanges'].items():
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound

    # Set atmospheric conditions
    reaction_oxygen = model.reactions.get_by_id('EX_o2_e')
    if atmosphere_type == 'aerobic':
        reaction_oxygen.lower_bound = -20
    elif atmosphere_type == 'anaerobic':
        reaction_oxygen.lower_bound = 0
    elif atmosphere_type is not None:
        raise ValueError

    return model
