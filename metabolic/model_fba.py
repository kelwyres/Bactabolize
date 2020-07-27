import gzip
import pathlib
import json
import re
import sys
import warnings


import cobra.io


from . import media_definitions


def run(model_fp, fba_types, spec_fp):
    print('\n========================================')
    print('running FBA')
    print('========================================')
    # Read in model
    with model_fp.open('r') as fh:
        model = cobra.io.load_json_model(fh)
    # Run FBA
    if not fba_types or 'individual' in fba_types:
        fba_individal_sources(model)
    if not fba_types or 'media' in fba_types:
        fba_media(model)
    if spec_fp and (not fba_types or 'spec' in fba_types):
        fba_spec(model, spec_fp)


def fba_spec(model, spec_fp):
    with spec_fp.open('r') as fh:
        spec_reaction_bounds = json.load(fh)
    for spec, reaction_bounds in spec_reaction_bounds.items():
        run_fba(model, reaction_bounds)


def fba_individal_sources(model):
    # Get formulas of metabolites involved in extracellular exchange
    metabolite_ids = list()
    for reaction in model.exchanges:
        for m in reaction.metabolites:
            [mid] = m.annotation['metanetx.chemical']
            metabolite_ids.append(mid)
    metabolite_formulas = get_formulas(metabolite_ids)

    # Distinct reactions
    # NOTE: these were defined from the described methods in the Shigella paper
    reaction_sources_distinct = {
        'EX_pi_e': 'inorganic phosphate',
        'EX_so4_e': 'inorganic sulfate',
        'EX_nh4_e': 'ammonia',
    }

    # Set up regex for selecting source type
    # NOTE: checking for single element characters insufficient as element symbol domain is degenerative
    element_base_re = '^.*%s([0-9A-Z]+.*)?$'
    carbon_re = re.compile(element_base_re % 'C')
    phosphate_re = re.compile(element_base_re % 'P')
    nitrogen_re = re.compile(element_base_re % 'N')
    sulfur_re = re.compile(element_base_re % 'S')

    # Get reactions and categories for FBAs
    fba_data = list()
    for reaction in model.exchanges:
        # Copy model and get metabolite formulas
        mformulas = list()
        for m in reaction.metabolites:
            [mid] = m.annotation['metanetx.chemical']
            mformulas.extend(metabolite_formulas[mid])
        # Get list of categories and additional reactions to enable
        categories = list()
        reaction_list = [reaction.id]
        if reaction.id in reaction_sources_distinct:
            categories.append(reaction_sources_distinct[reaction.id])
            reaction_list.append('EX_glc__D_e')
        else:
            if any(carbon_re.match(formula) for formula in mformulas):
                categories.append('carbon')
            if any(phosphate_re.match(formula) for formula in mformulas):
                categories.append('phosphate')
            if any(nitrogen_re.match(formula) for formula in mformulas):
                categories.append('nitrogen')
            if any(sulfur_re.match(formula) for formula in mformulas):
                categories.append('sulfur')
            # Add the default carbon source if no other is present
            if 'carbon' not in categories:
                reaction_list.append('EX_glc__D_e')
        # If no reaction is not relevant, ignore
        if not categories:
            continue
        fba_data.append((reaction_list, categories))

    # TODO: record this information in some parsable format
    # Execute FBAs
    for reaction_list, categories in fba_data:
        # Set lower bounds for reactions
        reaction_bounds = {r: -1000 for r in reaction_list}
        run_fba(model, reaction_bounds)


def fba_media(model):
    # TODO: define other media and for each get (1) objective biomass (2) all fluxes
    # TODO: better to simulate growth on media by reducing rx bounds via model.medium?
    #           print(model.medium)
    #           medium = model.medium
    #           medium['EX_o2_e'] = 0.0
    #           model.medium = medium
    # TODO: enable iterating all defined media
    run_fba(model, media_definitions.m9)


def run_fba(model, reaction_bounds):
    for reaction in model.exchanges:
        reaction.lower_bound = 0
    for reaction_id, lower_bound in reaction_bounds.items():
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            msg = f'warning: model does not contain reaction {reaction_id}'
            print(msg, file=sys.stderr)
        reaction.lower_bound = lower_bound

    # Prevent warnings from cobrapy
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Attempt to optimise
        solution = model.optimize()
        if solution.status == 'infeasible':
            print('warning: infeasible solution', file=sys.stderr)
        else:
            print()
            print(model.summary())


def get_formulas(queries):
    # Set filepaths
    map_fp = pathlib.Path('data/compound_id_map.tsv')
    formula_fp = pathlib.Path('data/formula.tsv.gz')

    # Get updated identifiers
    id_map = parse_id_map(map_fp)
    identifiers = {q: get_latest_id(q, id_map) for q in queries}

    # Get inverse dict
    identifiers_inv = dict()
    for query, nids in identifiers.items():
        for nid in nids:
            if nid not in identifiers_inv:
                identifiers_inv[nid] = list()
            identifiers_inv[nid].append(query)

    # Get formula
    formulas = {q: list() for q in queries}
    with gzip.open(formula_fp, 'rt') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for identifier, formula in line_token_gen:
            if identifier in identifiers_inv:
                for query in identifiers_inv[identifier]:
                    formulas[query].append(formula)
    return formulas


def get_latest_id(identifier, id_map, *, version='4.0'):
    if identifier not in id_map:
        return [identifier]
    id_out = list()
    id_queue = [identifier]
    while id_queue:
        cid = id_queue.pop(0)
        data = id_map[cid]
        for nid, vid in data:
            if vid == version:
                id_out.append(nid)
            else:
                id_queue.append(nid)
    return(id_out)


def parse_id_map(filepath):
    id_map = dict()
    with filepath.open('r') as fh:
        line_token_gen = (line.rstrip().split('\t') for line in fh)
        for id1, id2, version in line_token_gen:
            if id1 not in id_map:
                id_map[id1] = list()
            id_map[id1].append((id2, version))
    return id_map
