import gzip
import pathlib
import re
import warnings
import sys


# Regex for selecting source type
element_base_re = '^.*%s([0-9A-Z]+.*)?$'
carbon_re = re.compile(element_base_re % 'C')
phosphate_re = re.compile(element_base_re % 'P')
nitrogen_re = re.compile(element_base_re % 'N')
sulfur_re = re.compile(element_base_re % 'S')


def prepare_element_source_data(model):
    fba_data = dict()
    for reaction in model.exchanges:
        # Set categories for metabolite
        assert len(reaction.metabolites) == 1
        [metabolite] = reaction.metabolites
        categories = list()
        if carbon_re.match(metabolite.formula):
            categories.append('carbon')
        if phosphate_re.match(metabolite.formula):
            categories.append('phosphate')
        if nitrogen_re.match(metabolite.formula):
            categories.append('nitrogen')
        if sulfur_re.match(metabolite.formula):
            categories.append('sulfur')
        # Skip metabolites that are not potential element sources
        if not categories:
            continue
        fba_data[reaction.id] = categories
    return fba_data

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


def run_fba(model, exchange_bounds):
    # Set reaction lower bounds
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
