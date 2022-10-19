<img src="https://user-images.githubusercontent.com/19924405/193505313-edd9453a-e4eb-4730-81b1-a2bd9e652721.png" width="50%">

A high-throughput genome-scale Bactabolize model construction pipeline. Bactabolize allows you to provide an input
genome (annotated or unannotated) and construct a strain-specific Bactabolize model using a reference model. Growth
experience such as Flux Balance Analysis (FBA) and single gene knockout analysis can then be performed on the models
under a variety of growth conditions and mediums.

## Table of contents

* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [Model construction](#model-construction)
* [Growth profiles and Flux Balance Analysis](#growth-profiles-and-flux-balance-analysis)
* [Troubleshooting models](#troubleshooting-models)
* [Reference models](#reference-models)
* [Requirements](#requirements)
* [Development](#requirements)
* [Citation](#citation)
* [License](#license)

## Quickstart

```bash
# Create and activate environment
conda create -n bactabolize -y -c scwatts -c bioconda -c conda-forge 'bactabolize==0.1.0'
conda activate bactabolize

# Create draft model
bactabolize draft_model \
    --assembly_fp input_assembly.fasta \
    --ref_genes_fp reference_model_genes.ffn \
    --ref_proteins_fp reference_model_genes.faa \
    --ref_model_fp reference_model.json \
    --biomass_reaction_id YOUR_BIOMASS_ID \
    --media_type m9 \
    --atmosphere_type aerobic \
    --min_coverage 25 \
    --min_pident 80 \
    --output_fp input_assembly_model

# Run FBA
bactabolize fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp FBA_spec_files/M9_media_spec.json \
    --output_fp input_assembly_model_FBA.tsv \
    --fba_open_value -20
```

## Model construction

Metabolic model construction is run using the `bactabolize draft_model` command. Once a model is constructed,
Bactabolize then tests the model for growth on your choice of media under your choice of atmosphere. If the model does not grow under these conditions, `bactabolize patch_model` should be run to add additional reactions.

### Options

#### Required

`--assembly_fp` - Input assembly for which a Bactabolize model will be generated. This can be either an unannotated
**fasta** file or an annotated **genbank** file. Bactabolize will honour the genbank annotations. **IMPORTANT:**
Recommended minimum assembly quality: ≤200 assembly graph dead ends (calculate from .gfa or fastg). If only contigs are
available, ≤130 contigs.

`--output_fp` - Output filename

`--ref_model_fp` - Reference of Bactabolize model in .json format.

`--biomass_reaction_id` - ID of reference model's Biomass function. DEFAULT: BIOMASS_

The reference genome data used to generate the input assembly model can be provided as nucleotide (.ffn) AND protein
multifasta (.faa) files. Useful if reference is a pan-model or multi-strain model and does not exist in a traditional
genbank format.

`--ref_genes_fp` - Reference genes in nucleotide sequence (fasta). Corresponds to `--ref_proteins_fp` and
`--ref_model_fp`

`--ref_proteins_fp` - Reference protein in amino acid sequence (fasta). Corresponds to `--ref_genes_fp` and
`--ref_model_fp`

Alternatively, reference genome data can be provided by a genbank file. Useful if reference is a single-strain model.

`ref_genbank_fp` - Reference genome (genbank). Corresponds to `--ref_model_fp`

#### Optional

`--media_type` - Choose growth media for model building. One of: cdm_mendoza, bg11, lb, lb_carveme, m9, nutrient, pmm5_mendoza, pmm7_mendoza, tsa,
tsa_sheep_blood. DEFAULT: m9

`--atmosphere_type` - Choose atmosphere for model building. One of: aerobic, anaerobic. DEFAULT: aerobic

`--min_coverage` - Set minimum query coverage percentage for bi-directional best hit for ortholog identification.
DEFAULT: 25

`--min_pident` - Set minimum identity percentage for bi-directional best hit for ortholog identification. DEFAULT: 80

`--min_ppos` - Set minimum protein similarity (positives) percentage for bi-directional best hit for ortholog
identification. DEFAULT: OFF. Can be used instead of `--min_pident` to allow for greater tolerance of
similarly-functional but different residues.

`--memote_report_fp` - MEMOTE model quality report output filepath. Note that this will add >5 minutes of compute time
PER assembly

`--no_reannotation` - Will prevent the re-annotation of the input genbank file with prodigal. DEFAULT: Off.

### Examples

```bash
# Create draft model for genbank input assembly using genbank reference on M9 media
# under aerobic conditions, at 25% query coverage and 80% protein similarity
bactabolize draft_model \
    --assembly_fp input_assembly.gbk \
    --ref_genbank_fp reference.gbk \
    --ref_model_fp reference_model.json \
    --biomass_reaction_id biomass_equation \
    --media_type m9 \
    --atmosphere_type aerobic \
    --output_fp input_assembly_model_qc_25_sim_85 \
    --min_coverage 25 \
    --min_ppos 80

# Create draft model for fasta input assembly using multifasta reference, anerobically on PMM7 media
# at 25% query coverage and 75% protein identity. Produce MEMOTE report
bactabolize draft_model \
    --assembly_fp input_assembly.fasta \
    --ref_genes_fp reference_model_genes.ffn \
    --ref_proteins_fp reference_model_genes.faa \
    --ref_model_fp reference_model.json \
    --biomass_reaction_id BIOMASS \
    --media_type pmm7 \
    --atmosphere_type anaerobic \
    --output_fp input_assembly_model_qc_25_sim_85 \
    --min_coverage 25 \
    --min_pident 75 \
    --memote_report_fp input_report
```

### Draft model outputs

| Filename                      | Description                                      |
| ---------                     |---------                                         |
| `assembly_id`.gbk             | Prodigal annotation of input assembly            |
| `assembly_id`\_model.json     | Metabolic model in .json format                  |
| `assembly_id`\_model.xml      | Metabolic model in .xml ([SMBL Level 3 Version 1)](https://co.mbine.org/specifications/sbml.level-3.version-1.core.release-1) |
| `assembly_id`\_model.html     | [MEMOTE](https://github.com/opencobra/memote) model report                             |

## Growth profiles and Flux Balance Analysis

Once a model has been generated, you can then test its growth profiles across a range of nutrient sources. This is
performed using Flux Balance Analysis (FBA) and will test each nutrient source under aerobic and anaerobic conditions.

FBA here has been designed around (1) simulating growth on various media, and (2) identifying extracellular
metabolites as a source of carbon, phosphorus, nitrogen, and/or sulfur. An extracellular metabolite is considered a potential
source of an element if it is present in the metabolite chemical formula (e.g. a metabolite is a potential carbon
source if it contains carbon). In order to test whether a metabolite can be utilised in a draft model as an element
source, FBA is performed with only the media-defined exchanges and the target metabolite exchange enabled. Additionally,
the exchange of any default element source is disable in the minimal media as required for testing (i.e. when testing a
metabolite as a source of carbon, the default carbon source exchange is disabled). Where a metabolite contains more than
one of the four elements, all combinations are tested.

### Options

#### Required

`--model_fp` - Input model (.json) to perform FBA on

`--fba_spec_fp` - FBA spec file. Example and explanation [here](#fba-spec-file).

`--output_fp` - Output filename for [FBA results](#fba-output-file-format) (tab-delimited)

#### Optional

`--fba_open_value` - Set objective value for nutrient sources tested during FBA. Should be a negative value
between -1 and -1000. -10 or -20 is probably most reasonable. DEFAULT: -1000

### Examples

```bash
# Produce FBA on input on M9 minimal media with an objective value of -20
bactabolize fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp data/fba_specs/m9_spec.json \
    --output_fp input_assembly_model_FBA.tsv \
    --fba_open_value -20

# Produce FBA on input on TSA media with an objective value of -10
bactabolize fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp data/fba_specs/tsa_spec.json \
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

The `*_spec.json` file defines the media, atmosphere and default element sources
FBA will use as a growth environment. Example:

```json
{
  "M9": {
    "fba_type": [
        "defined_exchanges_only",
        "potential_element_sources"
    ],
    "atmosphere": ["aerobic", "anaerobic"],
    "media_type": "m9",
    "default_element_sources": {
      "carbon": "EX_glc__D_e",
      "phosphorus": "EX_pi_e",
      "nitrogen": "EX_nh4_e",
      "sulfur": "EX_so4_e"
    }
  }
}
```

m9_spec.json details: The `fba_type` field sets the type of FBA to run; both
simple media assessment and identification of potential element sources will be done here. The `atmosphere` allows
aerobic and anerobic conditions to be tested. The `media_type` field specifies the growth media (m9). 
`default_element_sources` specifies the default exchanges to replace when assessing potential element sources.

See `data/fba_specs` directory to find commonly used bacterial medias, or add custom ones.


## Troubleshooting models

The ability of a draft model to produce biomass on minimal media is assessed during creation. When a model fails this
test, troubleshooting information describing metabolites, reactions, and genes required to produce biomass is written to
disk.

In order to fix the model you must then add these missing reactions into a 'patch' file. Only the reactions are
required, their metabolites and associated genes (if any) will be added automatically. Then the `bactabolize
patch_model` command can be run, which will perform targeted gap-filling to repair the model. Example `patch.json` file:

```json
{
  "K_variicola_variicola_342": {
    "reactions": {
      "TDPDRE": "add",
      "TDPDRR": "add"
    }
  }
}
```

This patch file specifies that the `K_variicola_variicola_342` model requires two reactions to be added.

### patch_model requirements

* Add missing `reactions` to patch file
* Make sure the model name, in this case, `K_variicola_variicola_342`, matches the model name
`"id":"K_variicola_variicola_342"`, found just above the `"compartments":{` line in the `.json` model file

```bash
# K_variicola_variicola_342.json fails to produce biomass as it lacks lacks DTDP-4-dehydrorhamnose 3,5-epimerase and
# DTDP-4-dehydrorhamnose reductase, which are required to create DTDP-L-rhamnose
# Here we assume DTDP-L-rhamnose is not essential for growth and patch the model accordingly
bactabolize patch_model \
    --draft_model_fp K_variicola_variicola_342.json \
    --ref_model_fp iYL1228_annotated.json \
    --patch_fp patch.json \
    --media_type m9 \
    --atmosphere_type aerobic \
    --output_fp K_variicola_variicola_342_patched.json

# Assess model with FBA
bactabolize fba \
    --model_fp K_variicola_variicola_342_patched.json \
    --fba_spec_fp data/fba_specs/m9_spec.json \
    --output_fp K_variicola_variicola_342_patched_FBA.tsv \
    --fba_open_value -20
```

#### Options

##### Required

`--draft_model_fp` - Input draft model which requires gap-filling (.json)

`--ref_model_fp` - Reference model to gap-fill the draft model (.json)

`--patch_fp` - Missing reactions to add to input draft model (.json)

`--output_fp` - Output filename

##### Optional

`--media_type` - Choose growth media for model building. One of: bg11, lb_carveme, lb, m9, nutrient, tsa,
tsa_sheep_blood. DEFAULT: m9

`--atmosphere_type` - Choose atmosphere for model building. One of: aerobic, anaerobic. DEFAULT: aerobic

`--memote_report_fp` - MEMOTE model quality report output filepath. Note that this will add >5 minutes of compute time
PER assembly

## Reference models

To run individual FBA on extracellular metabolites, they must be annotated with the respective chemical formula in the
model. If you have a model that contains `metanetx` identifiers for metabolites (i.e. a BiGG model), you can add
metabolite formulas using the [`BiGG model compound
annotator`](https://github.com/scwatts/bigg_model_compound_annotator).

## Requirements

* python ==3.9
* biopython 1.79
* blast 2.12.0
* cobra 0.21.0
* prodigal 2.6.3
* memote 0.13.0

## Citation

Please cite the Bactabolize and COBRApy papers if you make use of Bactabolize

* bioRxiv and later, published paper here
* Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC
  Syst Biol 7, 74 (2013). <https://doi.org/10.1186/1752-0509-7-74>

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
which -a bactabolize
bactabolize --version
```

## License

[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
