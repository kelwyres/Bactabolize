##### Adding additional annotation data to COBRA metabolic models

"""
This script improves metabolic models by adding additional annotations for model metabolites and reactions by using a reference model (universal model from BiGG - https://zenodo.org/record/6713279 or a user-input model included).
User can provide reference and input model files in any format (.xml, .json or .sbml).
"""

### Authors
"""
Jon M Monk (https://scholar.google.com/citations?user=wJ4CAlwAAAAJ&hl=en)
Ben Vezina (https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en)
"""

### Dependancies
"""
cobra=0.20.0
memote=0.13.0
urllib3=1.26.3
"""

### Import packages
import cobra
from cobra.io import read_sbml_model, load_model, load_json_model, save_json_model, write_sbml_model
from glob import glob
import argparse
import os
import shutil
import urllib3
from urllib.request import urlopen
import json
import memote

### Define user input arguments
parser = argparse.ArgumentParser(description = 'Help for improve_model_annotations.py')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-m', '--models', help = 'Directory of model files to be annotated. Models can be in .xml, .json or .sbml formats', required = True) # args.models
required.add_argument('-r', '--reference', help = 'Reference model file annotations will be copied from. Reference model provided (iWFL1372_iML1515_iYS1720.json). Can also use BiGG Universal model found here: https://zenodo.org/record/6713279', required = True) # args.reference
optional.add_argument('-y', '--ya_name', help = 'What is ya name?') # args.ya_name
args = parser.parse_args()
# args = parser.parse_args('-m models -r reference/iWFL1372_iML1515_iYS1720.json'.split())

### Define functions

# Function for converting .json or SBML files to .xml to allow MEMOTE comparisons for both pre- and post-annotation 
def convert_to_xml(MODEL_DIR):
    """
    Before improving a model, the files are converted to SBML v3.1 format (.xml). This is so a MEMOTE score can be generated. Takes input model files
    """
    print("Converting input files to SBML v3.1 (.xml)")
    for model_file in MODEL_DIR:
        # Convert sbml old format to .xml
        if model_file.endswith('.sbml'):
            model, errors = validate_sbml_model(model_file) # doesn't print annoying errors due to legacy format...
            write_sbml_model(model, str(args.models) + '/' + str(model) + '_convert.xml')
            print("Saved legacy SBML as sbml v3.1 (.xml)")
        # read json files
        elif model_file.endswith('.json'):
        # if model_file.suffix == '.json': # use this for string-based method
            model = load_json_model(model_file)
            # Save as .xml for MEMOTE reports.
            write_sbml_model(model, str(args.models) + '/' + str(model) + '_convert.xml')
            print("Saved .json as sbml v3.1 (.xml)")

# Function for adding annotations to reactions and metabolites from a reference model 
def annotate(MODEL_FILES, REFERENCE):
    """
    This function adds annotations from the input reference model file. Takes input reference and a directory containing xml files 
    """
    # Run through all model files via loop
    for model_file in MODEL_FILES:
        # Define lists
        missing_metabolites = []
        missing_reactions = []
        # Load model
        print('loading %s'%model_file)
        model = cobra.io.read_sbml_model(model_file)
        for r in model.reactions:
            if r.id in REFERENCE.reactions:
                r.annotation = REFERENCE.reactions.get_by_id(r.id).annotation
            else:
                missing_reactions.append(r.id)
        for m in model.metabolites:
            if m.id in REFERENCE.metabolites:
                m.annotation = REFERENCE.metabolites.get_by_id(m.id).annotation
                print('\t downloading formula for %s'%m.id)
                try:
                    res = urllib3.urlopen('http://bigg.ucsd.edu/api/v2/universal/metabolites/' + str(m.id[:-2])) # see http://bigg.ucsd.edu/data_access
                    res = json.loads(res.read())
                    m.formula = res['formulae'][0]
                except:
                    print('error for %s'%m.id)
                    pass
            else:
                missing_metabolites.append(m.id)
        # Print summary of missing reactions and metabolites
        print ("Reactions with missing annotation data:\n" + str(set(missing_reactions)))
        print ("Metabolites with missing annotation data:\n" + str(set(missing_metabolites)))
        write_sbml_model(model, str(args.models) + '/' + str(model) + '_new_annotations.xml') # Print as xml - REQUIRED FOR MEMOTE
        save_json_model(model, str(args.models) + '/' + str(model) + '_new_annotations.json') # Print as json 
        print("Saved new model files with annotations")

# Function for producing MEMOTE reports via python to evaluate the new annotation scores
def MEMOTE_report_card(XML_MODELS):
    for model_file in XML_MODELS:
        model = read_sbml_model(model_file)
        result = memote.test_model(model, results=True)
        report = memote.snapshot_report(result[1], config=None, html=True)
        with open(str(model_file[:-4]) + "_report.html", "w") as handle:
            handle.write(report)

### Run script
if args.ya_name is None:
    print("Hi stranger, lets improve some models")
else:
    print("Hi " + str(args.ya_name) + ", lets improve some models")

# Read models in any format (.xml, .json, .sbml)
files = glob(str(args.models) + '/*.json') + glob(str(args.models) + '/*.xml') + glob(str(args.models) + '/*.sbml') # files = glob('models/*.json') + glob('models/*.xml') + glob('models/*.sbml')

# Run conversion module to xml
convert_to_xml(files)

# Load only xml files including converted ones
files_xml = glob(str(args.models) + '/*.xml')

# Define reference model
if args.reference.endswith('.sbml'):
    universal, errors = validate_sbml_model(args.reference) # doesn't print annoying errors due to legacy format
elif args.reference.endswith('.json'):
    universal = load_json_model(args.reference)
elif args.reference.endswith('.xml'):
    universal = read_sbml_model(args.reference)

# Annotate model using reference (or universal model)
annotate(files_xml, universal)

# Load only xml files again but include newly annotated files
files_pre_post = glob(str(args.models) + '/*.xml')

# Generate MEMOTE scores
MEMOTE_report_card(files_pre_post)
