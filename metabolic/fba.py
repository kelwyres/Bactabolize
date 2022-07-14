import re
import warnings
import sys


# Regex for selecting source type
element_base_re = '^.*%s([0-9A-Z]+.*)?$'
carbon_re = re.compile(element_base_re % 'C')
phosphorus_re = re.compile(element_base_re % 'P')
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
        if phosphorus_re.match(metabolite.formula):
            categories.append('phosphorus')
        if nitrogen_re.match(metabolite.formula):
            categories.append('nitrogen')
        if sulfur_re.match(metabolite.formula):
            categories.append('sulfur')
        # Skip metabolites that are not potential element sources
        if not categories:
            continue
        fba_data[reaction.id] = categories
    return fba_data


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
        return solution.objective_value
