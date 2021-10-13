# metabolic
A high-throughput metabolic model construction pipeline.


## Table of contents
* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [FBA overview](#fba-overview)
* [Troubleshootinga models](#troubleshooting-models)
* [Requirements](#requirements)
* [Outputs](#outputs)
* [License](#license)


## Quickstart
```bash
# Clone the repository
git clone https://github.com/scwatts/metabolic_pipeline.git && cd metabolic_modelling/

# Install dependencies and activate conda environment
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file misc/conda_packages.txt
conda activate $(pwd -P)/conda_env
pip install memote==0.13.0 # This takes ages to install. Don't worry, it's not you

# Create draft model
mkdir -p output/kvt/
./metabolic-runner.py draft_model \
    --assembly_fp data/K_variicola_tropicalensis_CDC4241-71.gbk \
    --ref_genbank_fp data/K_pneumoniae_MGH78578.gbk \
    --ref_model_fp data/iYL1228_annotated.json \
    --output_fp output/kvt/K_variicola_tropicalensis_CDC4241-71.json

# Run FBA
./metabolic-runner.py fba \
    --model_fp output/kvt/K_variicola_tropicalensis_CDC4241-71_model.json \
    --fba_spec_fp data/m9_fba_spec.json \
    --output_fp output/kvt/K_variicola_tropicalensis_CDC4241-71_fba.tsv
```

## FBA overview
Flux balance analysis here has been designed around (1) simulating growth on various media, and (2) identifying extracellular
metabolites as a source of carbon, phosphate, nitrogen, and/or sulfur. An extracellular metabolite is considered a potential
source of an element if it is present in the metabolite chemical formular (e.g. a metabolite is a potential carbon source if
it contains carbon). In order to test whether a metabolite can be utilised in a draft model as an element source, FBA is
performed with only the media-defined exchanges and the target metabolite exchange enabled. Additionally, the exchange of
any default element source is disable in the minimal media as required for testing (i.e. when testing a metabolite as a
source of carbon, the default carbon source exchange is disabled). Where a metabolite contains more than one of the four
elements, all combinations are tested. The `fba_spec` JSON file defines the type of FBA to perform, the media that it is
performed on, and the default element sources. Here is an example:
```json
{
  "m9": {
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
      "EX_sel_e":      -1000,
      "EX_slnt_e":     -1000,
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
This spec defines a single media setting, `m9`, to perform FBA. The `fba_type` field sets the type of FBA to run; both simple
media assessment and identification of potential element sources will be done here. The `exchanges` field specifies the media
and `default_element_sources` specifies the default exchanges to use when assessing potential element sources.

Information on FBA output can be found [here](#fba-results-file-format).



## Troubleshooting models
The ability of a draft model to produce biomass on minimal media is assessed during creation. When a model fails this test,
troubleshooting information describing metabolites, reactions, and genes required to produce biomass is written to disk.

In order to fix the model you must first determine what changes must be made to the model and then transcribe those changes
into a 'patch' file. Here is an example patch file:
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
This patch file specifies that the `K_variicola_variicola_342` model requires two reactions to be added. The `patch` command
is used to apply these changes. A brief example is shown below:
```bash
# Generate a draft model
mkdir -p output/kvv/
./metabolic-runner.py draft_model \
    --assembly_fp data/K_variicola_variicola_342.gbk \
    --ref_genbank_fp data/K_pneumoniae_MGH78578.gbk \
    --ref_model_fp data/iYL1228_annotated.json \
    --output_fp output/kvv/K_variicola_variicola_342.json

# This model fails to produce biomass as it lacks lacks DTDP-4-dehydrorhamnose 3,5-epimerase and
# DTDP-4-dehydrorhamnose reductase, which are required to create DTDP-L-rhamnose
# Here we assume DTDP-L-rhamnose is not essential for growth and patch the model accordingly
./metabolic-runner.py patch_model \
    --draft_model_fp output/kvv/K_variicola_variicola_342_model.json \
    --ref_model_fp data/iYL1228_annotated.json \
    --patch_fp data/kvv_patch_biomass.json \
    --output_fp output/kvv/K_variicola_variicola_342_model_patched.json

# Assess model with FBA
./metabolic-runner.py fba \
    --model_fp output/kvv/K_ariicola_variicola_342_model_patched.json \
    --fba_spec_fp data/m9_fba_spec.json \
    --output_fp output/kvv/K_variicola_variicola_342_fba.tsv
```

## Requirements
### Reference model
To run individual FBA on extracellular metabolites, they must be annotated with the respective chemical formula in the model.
If you have a model that contains `metanetx` identifiers for metabolites (i.e. a BiGG model), you can add metabolite formulas
using the [`BiGG model compound annotator`](https://github.com/scwatts/bigg_model_compound_annotator).


## Outputs
| Filename                      | Description                           |
| ---------                     |---------                              |
| `assembly_id`.gbk             | Prodigal annotation of input assembly |
| `assembly_id`\_model.json     | Metabolic model                       |
| `assembly_id`\_fba.tsv        | Results of FBA                        |


### FBA results file format
| Column name       | Description                                                                   |
| ---------         |---------                                                                      |
| `fba_type`        | Name of the FBA run (defined exchanges only, or all potential element sources |
| `spec_name`       | FBA specification name                                                        |
| `atmosphere`      | Type of atmosphere (aerobic: O2; anaerobic: no O2)                            |
| `exchange`        | Name of exchange assessed                                                     |
| `categories`      | Elements for which the exchange was assessed as being a potential source      |
| `objective value` | Biomass objective value                                                       |


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
