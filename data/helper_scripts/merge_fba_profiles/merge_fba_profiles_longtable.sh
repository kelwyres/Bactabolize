#!/bin/bash

### ------------------About the script-------------------

# This script will merge multiple Flux Balance Analysis (FBA) profiles which are output from Bactabolize's fba module. Input files are *.tsv all within the same directory

### ---------------------Usage-----------------------

# bash merge_fba_profiles_longtable.sh -f fba_profiles_directory

### ------------------Dependencies------------------

# None

### ------------------Installation------------------

# None required. Run tool via bash command line. 

### --------------Create user options---------------

# Define help section
COMMAND_LINE_OPTIONS_HELP="
This script will merge multiple Flux Balance Analysis (FBA) profiles which are output from Bactabolize's fba module. Input files are *.tsv all within the same directory

Usage: bash merge_fba_profiles_longtable.sh [-h] -f fba_profiles_directory

Command line options:
    Required:
        -f          Directory containing FBA profiles outputted from Bactabolize's fba command. Files in .tsv format.
		
    Optional:
        -h          Print this help menu

Example:
	bash merge_fba_profiles_longtable.sh -f fba_profiles_directory
"

# Define options
VALID_COMMAND_LINE_OPTIONS="f:h"

while getopts $VALID_COMMAND_LINE_OPTIONS options; do
    #echo "option is " $options
    case $options in
        h)
            echo "$COMMAND_LINE_OPTIONS_HELP"
            exit $E_OPTERROR;
        ;;
        f)
            fba_profiles_directory=${OPTARG};;
        \?)
            echo "bash merge_fba_profiles_longtable.sh -h for help";
            echo "$COMMAND_LINE_OPTIONS_HELP"
            exit $E_OPTERROR;
        ;;
    esac
done

# Test to see if the model flag is empty
if [ -z "$fba_profiles_directory" ]
then
   echo "$COMMAND_LINE_OPTIONS_HELP"
fi

### -------------------Run script--------------------

# For loop processing each FBA profile from Bactabolize by adding isolate as col 1
for tsv in $fba_profiles_directory/*.tsv; do
tag="$(basename $tsv .tsv)" # basename removes path and you can specify the tag
# this line adds a new col 1 with the isolate name, then adds 'Isolate' as a header.
sed -e "s/^/$tag\t/g" $tsv | sed -e "1s/$tag/"Isolate"/" > $fba_profiles_directory/"$tag"_alt.tsv
done

# Merge all together and use only header from first file! https://apple.stackexchange.com/questions/80611/merging-multiple-csv-files-without-merging-the-header
awk '(NR == 1) || (FNR > 1)' $fba_profiles_directory/*_alt.tsv > $fba_profiles_directory/fba_profiles_merged_longtable.tsv
rm -rf $fba_profiles_directory/*_alt.tsv
