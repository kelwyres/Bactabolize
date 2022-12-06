# Improving metabolic model annotations

This script will improve the annotations and metadata tagging of metabolic model reactions and metabolites using a reference model of high quality. For most cases, this is the [universal.json](http://bigg.ucsd.edu/static/namespace/universal_model.json). For an updated universal.json compatible with modern releases of COBRApy, see [here](https://zenodo.org/record/6713279)

## Overview and process
1. The script takes input directory contaning metabolic model file/s and reference model file (model file formats: `.xml`, `.json`, `.sbml`) 
2. Convert model file to the `.xml` SBML v3.1 format
3. Use the reference model and the univeral model to improve the annotations of the input model
4. MEMOTE reports will be generated for the pre- and post-annotation model files for comparison.

## Quick start
```
python improve_model_annotations.py --models model_directory --reference reference/iWFL1372_iML1515_iYS1720.json
```

## Dependancies
- cobra=0.20.0
- memote=0.13.0
- urllib3=1.26.3
- python=3.6

## Installation
The best way to install this is to download this directory (including the `reference` directory) and run the script. Recommend to install within a conda environment.

```
# Create conda environment
conda create -y --name metabolic_model_env python=3.6.12

# Activate environment
conda activate metabolic_model_env

# Install dependancies
conda install -y -c bioconda cobra=0.20.0
conda install -y -c conda-forge urllib3
pip install memote==0.13.0

# Test script
python improve_model_annotations.py -h
```

## Usage
```
improve_model_annotations.py [-h] -m MODELS -r REFERENCE [-y YA_NAME]

required arguments:
  -m MODELS, --models MODELS
                        Directory of model files to be annotated
  -r REFERENCE, --reference REFERENCE
                        Reference model file annotations will be copied from. Reference model provided (iWFL1372_iML1515_iYS1720.json). Can also use BiGG Universal model found
                        here: https://zenodo.org/record/6713279

optional arguments:
  -y YA_NAME, --ya_name YA_NAME
                        What is ya name?
```

## Tips
The better the reference model you use, the better your model annotations will likely be. The iWFL1372_iML1515_iYS1720.json model is a [merged model](https://github.com/bananabenana/Metabolic_modelling_scripts/tree/main/merge_metabolic_models) of [iWFL1372](http://bigg.ucsd.edu/static/models/iWFL_1372.xml), [iML1515](http://bigg.ucsd.edu/static/models/iML1515.xml) and [iYS1720](http://bigg.ucsd.edu/static/models/iYS1720.xml), all of which are well-curated and have good model annotations. The universal.json model is another great choice.

## Authors

- Jon M Monk (https://scholar.google.com/citations?user=wJ4CAlwAAAAJ&hl=en)
- Ben Vezina (https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en)
