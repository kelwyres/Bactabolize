import gzip
import pathlib
import re
import warnings
import sys


# Default element sources
source_defaults = {
    'carbon': 'EX_glc__D_e',
    'phosphate': 'EX_pi_e',
    'sulfur': 'EX_so4_e',
    'nitrogen': 'EX_nh4_e'
}

# Regex for selecting source type
element_base_re = '^.*%s([0-9A-Z]+.*)?$'
carbon_re = re.compile(element_base_re % 'C')
phosphate_re = re.compile(element_base_re % 'P')
nitrogen_re = re.compile(element_base_re % 'N')
sulfur_re = re.compile(element_base_re % 'S')


def get_individual_data(model, metabolite_formulas):
    # Get reactions and categories for FBAs
    # Additionally record list of exchanges to block during FBA
    fba_data = list()
    exchange_sources = set()
    for reaction in model.exchanges:
        # Get metabolite formulas
        mformulas = list()
        for m in reaction.metabolites:
            [mid] = m.annotation['metanetx.chemical']
            mformulas.extend(metabolite_formulas[mid])
        # Categories
        categories = list()
        if any(carbon_re.match(formula) for formula in mformulas):
            categories.append('carbon')
        if any(phosphate_re.match(formula) for formula in mformulas):
            categories.append('phosphate')
        if any(nitrogen_re.match(formula) for formula in mformulas):
            categories.append('nitrogen')
        if any(sulfur_re.match(formula) for formula in mformulas):
            categories.append('sulfur')
        # If exchange cannot possibly act as a source for one of the four elements, skip
        if not categories:
            continue
        # Reactions
        reaction_list = [reaction.id]
        source_missing = [s for s in source_defaults if s not in categories]
        for source in source_missing:
            reaction_list.append(source_defaults[source])
        exchange_sources.add(reaction.id)
        # Add both aerobic and anaerobic conditions
        for atmo, reaction in (('anaerobic', 'EX_co2_e'), ('aerobic', 'EX_o2_e')):
            categories_atmo = [*categories, atmo]
            reaction_list_atmo = [*reaction_list, reaction]
            fba_data.append((reaction_list_atmo, categories_atmo))
    # Add CO2 and O2 exchanges to block list
    exchange_sources.add('EX_co2_e')
    exchange_sources.add('EX_o2_e')
    return fba_data, exchange_sources


def run_fba(model, exchange_bounds, *, exchange_block=None):
    # Set reaction lower bounds
    block_list = exchange_block if exchange_block else model.exchanges
    for reaction in block_list:
        if isinstance(reaction, str):
            reaction = model.reactions.get_by_id(reaction)
        reaction.lower_bound = 0
    for reaction_id, lower_bound in exchange_bounds.items():
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
            return 'infeasible'
        else:
            return solution.objective_value


def get_formulas(queries):
    # Set filepaths
    base_dir = pathlib.Path(__file__).parent / 'data'
    map_fp = base_dir / 'compound_id_map.tsv'
    formula_fp = base_dir / 'formula.tsv.gz'
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
