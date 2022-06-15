# metabolic

A high-throughput genome-scale metabolic model construction pipeline. metabolic allows you to provide an input genome (annotated or unannotated) and construct a strain-specific metabolic model using a reference model. Growth experience such as Flux Balance Analysis (FBA) can then be performed on the models under a variety of growth conditions and mediums.

## Table of contents

* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [Model construction](#model-construction)
* [Growth profiles and Flux Balance Analysis](#growth-profiles-and-flux-balance-analysis)
* [Troubleshooting models](#troubleshooting-models)
* [Requirements](#requirements)
* [License](#license)

## Quickstart

```bash
### Installation

# Create environment
conda create -n -y metabolic_v_0.0.1

# Activate environment
conda activate metabolic_v_0.0.1

# Install
conda install -y -c scwatts -c bioconda -c conda-forge metabolic

### Construct and test metabolic models

# Create draft model
metabolic draft_model \
    --assembly_fp input_assembly.fasta \
    --ref_genes_fp reference_model_genes.ffn \
    --ref_proteins_fp reference_model_genes.faa \
    --ref_model_fp reference_model.json \
    --output_fp input_assembly_model \
    --min_coverage 25 \
    --min_pident 80

# Run FBA
metabolic fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp FBA_spec_files/M9_media_spec.json \
    --output_fp input_assembly_model_FBA.tsv \
    --fba_open_value -20
```

## Model construction

Metabolic model construction is run using the `metabolic draft_model` command. Once a model is constructed, metabolic then tests the model for growth on M9 minimal media with glucose. If the model does not grow under these conditions, `metabolic patch_model` should be run to add additional reactions.

### Options

`--assembly_fp` - (REQUIRED) Input assembly for which a metabolic model will be generated. This can be either an unannotated **fasta** file or an annotated **genbank** file. metabolic will honour the genbank annotations. Useful if you already have an annotation you want to generate a model from. 

`--output_fp` - (REQUIRED) Output filename

`--ref_model_fp` - (REQUIRED) Reference of metabolic model in .json format.

#### The reference genome data used to generate the input assembly model can be provided as nucleotide (.ffn) AND protein multifasta (.faa) files. Useful if reference is a pan-model or multi-strain model and does not exist in a traditional genbank format.

`--ref_genes_fp` - (REQUIRED) Reference genes in nucleotide sequence (fasta). Corresponds to `--ref_proteins_fp` and `--ref_model_fp`

`--ref_proteins_fp` - (REQUIRED) Reference protein in amino acid sequence (fasta). Corresponds to `--ref_genes_fp` and `--ref_model_fp` 


#### Alternatively, reference genome data can be provided by a genbank file. Useful if reference is a single-strain model.

`ref_genbank_fp` - (REQUIRED) Reference genome (genbank). Corresponds to `--ref_model_fp` 

#### Other options
`--min_coverage` - (OPTIONAL) Set minimum query coverage percentage for bi-directional best hit for ortholog identification. DEFAULT: 25

`--min_pident` (OPTIONAL) Set minimum identity percentage for bi-directional best hit for ortholog identification. DEFAULT: 80

`--min_ppos` (OPTIONAL) Set minimum protein similarity (positives) percentage for bi-directional best hit for ortholog identification. DEFAULT: OFF. Can be used instead of `--min_pident` to allow for greater tolerance of similarly-functional but different residues.

`--no_reannotation` - (OPTIONAL) Will prevent the re-annotation of the input genbank file with prodigal. DEFAULT: Off.

### Examples
```
# Create draft model for genbank input assembly using genbank reference, at 25% query coverage and 80% protein similarity 
metabolic draft_model \
    --assembly_fp input_assembly.gbk \
    --ref_genbank_fp reference.gbk \
    --ref_model_fp reference_model.json \
    --output_fp input_assembly_model_qc_25_sim_85 \
    --min_coverage 25 \
    --min_ppos 80
    
# Create draft model for fasta input assembly using multifasta reference, at 25% query coverage and 75% protein identity 
metabolic draft_model \
    --assembly_fp input_assembly.fasta \
    --ref_genes_fp reference_model_genes.ffn \
    --ref_proteins_fp reference_model_genes.faa \
    --ref_model_fp reference_model.json \
    --output_fp input_assembly_model_qc_25_sim_85 \
    --min_coverage 25 \
    --min_pident 75
```

### Draft model outputs

| Filename                      | Description                                      |
| ---------                     |---------                                         |
| `assembly_id`.gbk             | Prodigal annotation of input assembly            |
| `assembly_id`\_model.json     | Metabolic model in .json format                  |
| `assembly_id`\_model.xml      | Metabolic model in .xml ([SMBL Level 3 Version 1)](https://co.mbine.org/specifications/sbml.level-3.version-1.core.release-1) |
| `assembly_id`\_model.html     | [MEMOTE](https://github.com/opencobra/memote) model report                             |


## Growth profiles and Flux Balance Analysis

Once a model has been generated, you can then test its growth profiles across a range of nutrient sources. This is performed using Flux Balance Analysis (FBA) and will test each nutrient source under aerobic and anaerobic conditions.

FBA here has been designed around (1) simulating growth on various media, and (2) identifying extracellular
metabolites as a source of carbon, phosphate, nitrogen, and/or sulfur. An extracellular metabolite is considered a potential
source of an element if it is present in the metabolite chemical formula (e.g. a metabolite is a potential carbon
source if it contains carbon). In order to test whether a metabolite can be utilised in a draft model as an element
source, FBA is performed with only the media-defined exchanges and the target metabolite exchange enabled. Additionally,
the exchange of any default element source is disable in the minimal media as required for testing (i.e. when testing a
metabolite as a source of carbon, the default carbon source exchange is disabled). Where a metabolite contains more than
one of the four elements, all combinations are tested.


### Options

`--model_fp` - (REQUIRED) Input model (.json) to perform FBA on

`--fba_spec_fp` - (REQUIRED) FBA spec file. Example and explanation [here](#fba-spec-file).

`--output_fp` - (REQUIRED) Output filename for [FBA results](#fba-output-file-format) (tab-delimited)

`--fba_open_value` - (OPTIONAL) Set objective value for nutrient sources tested during FBA. Should be a negative value between -1 and -1000. -10 or -20 is probably most reasonable. DEFAULT: -1000

### Examples

```
# Produce FBA on input on M9 minimal media with an objective value of -20
metabolic fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp FBA_spec_files/M9_media_spec.json \
    --output_fp input_assembly_model_FBA.tsv \
    --fba_open_value -20
    
# Produce FBA on input on TSA media with an objective value of -10
metabolic fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp FBA_spec_files/TSA_media_spec.json \
    --output_fp input_assembly_model_FBA.tsv \
    --fba_open_value -10 
```

### FBA output file format
A `assembly_id`\_fba.tsv tab-delimited file will be produced:

| Column name       | Description                                                                   |
| ---------         |---------                                                                      |
| `fba_type`        | Name of the FBA run (defined exchanges only, or all potential element sources |
| `spec_name`       | FBA specification name                                                        |
| `atmosphere`      | Type of atmosphere (aerobic: O2; anaerobic: no O2)                            |
| `exchange`        | Name of exchange assessed                                                     |
| `categories`      | Elements for which the exchange was assessed as being a potential source      |
| `objective value` | Biomass objective value                                                       |


### FBA spec file
The `*_spec.json` file defines the type of FBA to perform, the
media that it is performed on, and the default element sources. Here is an example:

```json
{
  "M9": {
    "fba_type": [
        "defined_exchanges_only",
        "potential_element_sources"
    ],
    "exchanges": {
      "EX_ca2_e":      -1000,
      "EX_cbl1_e":     -0.01,
      "EX_cl_e":       -1000,
      "EX_cobalt2_e":  -1000,
      "EX_cu2_e":      -1000,
      "EX_fe2_e":      -1000,
      "EX_fe3_e":      -1000,
      "EX_glc__D_e":   -20,
      "EX_h2o_e":      -1000,
      "EX_h_e":        -1000,
      "EX_k_e":        -1000,
      "EX_mg2_e":      -1000,
      "EX_mn2_e":      -1000,
      "EX_mobd_e":     -1000,
      "EX_na1_e":      -1000,
      "EX_nh4_e":      -1000,
      "EX_ni2_e":      -1000,
      "EX_pi_e":       -1000,
      "EX_so4_e":      -1000,
      "EX_tungs_e":    -1000,
      "EX_zn2_e":      -1000
    },
    "default_element_sources": {
      "carbon": "EX_glc__D_e",
      "phosphate": "EX_pi_e",
      "nitrogen": "EX_nh4_e",
      "sulfur": "EX_so4_e"
    }
  }
}
```
This spec defines a single media setting, `M9`, to perform FBA. The `fba_type` field sets the type of FBA to run; both simple
media assessment and identification of potential element sources will be done here. The `exchanges` field specifies the media
and `default_element_sources` specifies the default exchanges to use when assessing potential element sources. 7 common bacterial medias (including TSA, LB, nutrient media, BG11 etc) have been included in the `FBA_spec_files` directory.

## Troubleshooting models

The ability of a draft model to produce biomass on minimal media is assessed during creation. When a model fails this test,
troubleshooting information describing metabolites, reactions, and genes required to produce biomass is written to disk.

In order to fix the model you must first determine what changes must be made to the model and then transcribe those changes
into a 'patch' file. Then the `metabolic patch_model` command can be run, which will perform targeted gap-filling to repair the model. 
Example `patch.json` file:

```json
{
  "K_variicola_variicola_342": {
    "reactions": {
      "TDPDRE": "add",
      "TDPDRR": "add"
    },
    "biomass_metabolites": {
    }
  }
}
```

This patch file specifies that the `K_variicola_variicola_342` model requires two reactions to be added.

```bash
# K_variicola_variicola_342.json fails to produce biomass as it lacks lacks DTDP-4-dehydrorhamnose 3,5-epimerase and
# DTDP-4-dehydrorhamnose reductase, which are required to create DTDP-L-rhamnose
# Here we assume DTDP-L-rhamnose is not essential for growth and patch the model accordingly
metabolic patch_model \
    --draft_model_fp K_variicola_variicola_342.json \
    --ref_model_fp iYL1228_annotated.json \
    --patch_fp patch.json \
    --output_fp K_variicola_variicola_342_patched.json

# Assess model with FBA
metabolic fba \
    --model_fp K_variicola_variicola_342_patched.json \
    --fba_spec_fp FBA_spec_files/M9_media.json \
    --output_fp K_variicola_variicola_342_patched_FBA.tsv \
    --fba_open_value -20
```

### Options

`--draft_model_fp` - (REQUIRED) Input draft model which requires gap-filling (.json)

`--ref_model_fp` - (REQUIRED) Reference model to gap-fill the draft model (.json)

`--patch_fp` - (REQUIRED) Missing reactions to add to input draft model (.json)

`--output_fp` - (REQUIRED) Output filename

### Reference models

To run individual FBA on extracellular metabolites, they must be annotated with the respective chemical formula in the model.
If you have a model that contains `metanetx` identifiers for metabolites (i.e. a BiGG model), you can add metabolite formulas
using the [`BiGG model compound annotator`](https://github.com/scwatts/bigg_model_compound_annotator).


## Requirements

## Development

Set up development environment and install pre-commit hooks

```bash
conda env create -p $(pwd -P)/conda_env/ -f requirements-dev.yaml -y
conda activate ./conda_env/
pre-commit install
```

Install as editable python package

```bash
# Install
pip install -e .

# Check
which -a metabolic
metabolic --version
```

## License

[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
