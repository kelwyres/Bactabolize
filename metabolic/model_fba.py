import sys


import cobra.io


from . import media_definitions


def run(model_fp):
    print('\n========================================')
    print('running FBA')
    print('========================================')
    #import pickle
    with model_fp.open('r') as fh:
        model = cobra.io.load_json_model(fh)
    #with open(f'pickled/model_{model_fp.stem}.bin', 'wb') as fh:
    #    pickle.dump(model, fh)
    #with open(f'pickled/model_{model_fp.stem}.bin', 'rb') as fh:
    #    model = pickle.load(fh)

    # TODO: define other media and for each get (1) objective biomass (2) all fluxes

    # TODO: better to simulate growth on media by reducing rx bounds via model.medium?
    #print(model.medium)
    #medium = model.medium
    #medium['EX_o2_e'] = 0.0
    #model.medium = medium

    for reaction in model.reactions:
        if reaction.id.startswith('EX_'):
            reaction.lower_bound = 0
    for reaction_id, lower_bound in media_definitions.m9.items():
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: model {model_fp.stem} does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound

    solution = model.optimize()
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        print()
        print(model.summary())
