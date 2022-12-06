# Merge multiple FBA profiles into longtable format

This script takes a directory containing  .tsv files from the [Bactabolize](https://github.com/kelwyres/Bactabolize) [fba](https://github.com/kelwyres/Bactabolize/wiki/4.-Predict-growth-profiles) command and merges them together into a longtable format for analysis in a program like R. 

## Overview and process
1. The script takes input directory of growth profile files (file format: `.tsv`) 
2. Merges together with column header
3. Outputs new file `fba_profiles_merged_longtable.tsv`


## Quick start
```
bash merge_fba_profiles_longtable.sh -f fba_profiles_directory
```

## Dependancies
- unix environment with basename, dirname and getopts

## Installation
The best way to install this is to download this directory and run the script directly in command line.

## Usage
```
# bash merge_fba_profiles_longtable.sh [-h] -f fba_profiles_directory

Command line options:
    Required:
        -f          Directory containing FBA profiles outputted from Bactabolize's fba command. Files in .tsv format.
		
    Optional:
        -h          Print this help menu

Example:
	bash merge_fba_profiles_longtable.sh -f fba_profiles_directory
 ```


## Authors

- Ben Vezina (https://scholar.google.com/citations?user=Rf9oh94AAAAJ&hl=en)
