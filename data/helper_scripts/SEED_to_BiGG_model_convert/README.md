# Convert SEED ID to BiGG ID in metabolic model for downstream analysis

This script will do a string replace from "cpd00007" (SEED ID) to "o2" (BiGG ID) for oxygen so models can be initialised and FBA, single gene knockout analysis etc performed. ModelSEED and KBase models will not effectively run in Bactabolize. This is because it utilises the BiGG ID for oxygen (o2_c, o2_p and o2_e).

## Overview and process
1. The script takes input metabolic model file/s (model file formats: `.xml`, `.json`) 
2. String replace for SEED IDs incompatible with Bactabolize
3. Outputs new file which can be used with Bactabolize's main functions


## Quick start
```
bash SEED_to_BiGG_model_convert.sh -m input.xml
```

## Dependancies
- unix environment with basename, dirname and getopts

## Installation
The best way to install this is to download this directory and run the script directly in command line.

```
# Test script
bash SEED_to_BiGG_model_convert.sh -h
```

## Usage
```
bash SEED_to_BiGG_model_convert.sh [-h] -m model.xml

Command line options:
    Required:
        -m          Model file to be converted. Takes .xml or .json
		
    Optional:
        -h          Print this help menu

Example:
	bash SEED_to_BiGG_model_convert.sh -m model.xml
	bash SEED_to_BiGG_model_convert.sh -m model.json
 ```


## Authors

- Ben Vezina (https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en)
