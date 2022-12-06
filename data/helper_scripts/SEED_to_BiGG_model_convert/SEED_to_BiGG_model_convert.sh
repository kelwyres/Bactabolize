#!/bin/bash

### ------------------About the script-------------------

# ModelSEED and KBase models will not effectively run in Bactabolize. This is because it utilises the BiGG ID for oxygen (o2_c, o2_p and o2_e). This script will do a string replace from cpd00007 (SEED ID for oxygen) so models can be initialised and FBA, single gene knockout analysis etc performed.

### ---------------------Usage-----------------------

# bash SEED_to_BiGG_model_convert.sh -m iKp1289.xml 


### ------------------Installation------------------

# None required. Run on unix environment with basename, dirname and getopts

### -------------------Run script--------------------

# Define help section
COMMAND_LINE_OPTIONS_HELP='
This script will do a string replace from "cpd00007" (SEED ID) to "o2" (BiGG ID) for oxygen so models can be initialised and FBA, single gene knockout analysis etc performed.

Usage: bash SEED_to_BiGG_model_convert.sh [-h] -m model.xml

Command line options:
    Required:
        -m          Model file to be converted. Takes .xml or .json
		
    Optional:
        -h          Print this help menu

Example:
	bash SEED_to_BiGG_model_convert.sh -m model.xml
	bash SEED_to_BiGG_model_convert.sh -m model.json
'

# Define options
VALID_COMMAND_LINE_OPTIONS="m:h"

while getopts $VALID_COMMAND_LINE_OPTIONS options; do
    #echo "option is " $options
    case $options in
        h)
            echo "$COMMAND_LINE_OPTIONS_HELP"
            exit $E_OPTERROR;
        ;;
        m)
            model=${OPTARG};;
        \?)
            echo "Usage: bash SEED_to_BiGG_model_convert.sh -h for help";
            echo "$COMMAND_LINE_OPTIONS_HELP"
            exit $E_OPTERROR;
        ;;
    esac
done

# Test to see if the model flag is empty
if [ -z "$model" ]
then
   echo "$COMMAND_LINE_OPTIONS_HELP"
fi

# Create filename tag
tag=$(basename "${model%.*}")
dir=$(dirname "${model}")
ext=${model##*\.} 

# Run either on json or xml
if [[ $ext == xml ]]; then
	echo "$tag is xml. Fixing model now"
	sed 's/cpd00007_c0/o2_c/g' $model > "$dir"/"$tag"_o2_fix.xml
	echo "..."
	sed -i 's/cpd00007_e0/o2_e/g' "$dir"/"$tag"_o2_fix.xml
	echo "..."
	sed -i 's/cpd00007_p0/o2_p/g' "$dir"/"$tag"_o2_fix.xml
	echo "Done!"

elif [[ $ext == json ]]; then
	echo "$tag is json. Fixing model now"
	echo "..."
	sed 's/cpd00007_c0/o2_c/g' $model > "$dir"/"$tag"_o2_fix.json
	echo "..."
	sed -i 's/cpd00007_e0/o2_e/g' "$dir"/"$tag"_o2_fix.json
	echo "..."
	sed -i 's/cpd00007_p0/o2_p/g' "$dir"/"$tag"_o2_fix.json
	echo "Done!"

else
	echo "Model format not supported. Try converting with model_format_converter.py"
	
fi
