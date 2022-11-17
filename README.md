<img src="https://user-images.githubusercontent.com/19924405/193505313-edd9453a-e4eb-4730-81b1-a2bd9e652721.png" width="50%">

A high-throughput genome-scale metabolic reconstruction and growth simulation pipeline.

Bactabolize is designed for rapid generation of strain-specific metabolic reconstructions from bacterial genome data using the approach described in [Norsigian et al. Nature Protocols 2020](https://www.nature.com/articles/s41596-019-0254-3). It leverages the [COBRApy toolkit](https://opencobra.github.io/cobrapy/) and takes an input genome assembly (annotated or unannotated) to construct a strain-specific draft model by comparison to a reference (ideally a multi-strain or 'pan'-reference model). It also allows high-throughput growth phenotype simulation via Flux Balance Analysis (FBA) e.g. to predict substrate usage profiles and Single Gene Knockout analysis (SGK) to predict the impacts of single gene knockout mutations. These can be performed under a variety of growth conditions and mediums.


## Table of contents

* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [Model construction](#model-construction)
* [Growth profiles and Flux Balance Analysis](#growth-profiles-and-flux-balance-analysis)
* [Troubleshooting models](#troubleshooting-models)
* [Metabolite IDs](#metabolite-ids)
* [Requirements](#requirements)
* [Citation](#citation)
* [Development](#requirements)
* [License](#license)

## Quickstart

The easiest way to install Bactabolize is via a conda environment. You can install conda as descrbed [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Then follow the steps below to create a conda environment for Bactabolize. Alternatively, [mamba](https://anaconda.org/conda-forge/mamba) can be used in place of conda, which allows for a faster install.

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
    --output_fp input_assembly

# Run FBA
bactabolize fba \
    --model_fp input_assembly_model.json \
    --fba_spec_fp FBA_spec_files/M9_media_spec.json \
    --output_fp input_assembly_model_FBA.tsv \
    --fba_open_value -20
```

## Model construction

Metabolic reconstructions ('models') are generated via the `bactabolize draft_model` command. Once a model is constructed, Bactabolize will test the model for the ability to simulate growth on your choice of media under your choice of atmosphere (we recommend testing for conditions in which all of your strains are expected to grow e.g. for *Klebsiella pneumoniae* we use m9 minimal media plus aerobic atmosphere). If the model does not simulate growth under these conditions, Bactabolize will attempt to identify the essential missing reactions via the [COBRApy gap-filling](https://cobrapy.readthedocs.io/en/latest/gapfilling.html) function. `bactabolize patch_model` can be used to add these reactions into the model.

### Reference model choice

To accurately capture the metabolism of your intended isolates and build high-quality models, choosing a good reference model is essential. Many approaches will use single strain models, which may not have enough diversity to capture your intended isolate's reactions. Other approaches will use universal models, which can overestimate reaction presence in a genome. We have constructed a high-quality, curated [Klebsiella pneumoniae Species Complex (*Kp*SC) pan metabolic model](https://github.com/kelwyres/KpSC-pan-metabolic-model), built from 37 *Kp*SC genomes, for use as a *Klebsiella* reference model.

### Growth mediums

You can choose your growth medium in which to simulate growth of your model. Bactabolize comes with a list of pre-generated medias in the `data/fba_specs/` directory, but users can make custom medias as well.

| Media type | Reference | File with ingredients |
|-------|-----------|-----------------------|
|   cdm_mendoza    |  [Mendoza *et al*., (2019)](https://doi.org/10.1186/s13059-019-1769-1)         |    cdm_mendoza_spec.json                   |
|    bg11   | This study, [ThermoFisher](https://www.thermofisher.com/au/en/home/technical-resources/media-formulation.353.html)           |  bg11_spec.json  |
|   lb    | This study. Peptone: [BD, 2015](https://legacy.bd.com/ds/technicalCenter/misc/lcn01558-bionutrients-manual.pdf), [Loginova *et al*., (1974)](https://doi.org/10.1007/BF00777001), [ThermoFisher, (2019)](https://assets.thermofisher.com/TFS-Assets/BPD/brochures/peptones-supplements-feeds-technical-reference-guide.pdf). Yeast extract: [Tomé, (2021)](https://doi.org/10.1021/acsfoodscitech.0c00131), [Plata *et al*., (2013)](https://doi.org/10.1007/s00216-013-7239-9), [](), [Liu *et al*., (2018)](https://doi.org/10.1016/j.ijbiomac.2018.06.145), [Blagović *et al*., (2001)](https://www.bib.irb.hr/94935?&rad=94935), [Blagović *et al*., (2005)](https://doi.org/10.1007/BF02931290), [Avramia *et al*., (2021)](https://doi.org/10.3390/ijms22020825)            |          lb_spec.json             |
|  lb_carveme     |  [Machado *et al*., (2018)](https://doi.org/10.1093/nar/gky537)         |          lb_carveme_spec.json             |
| m9      |   [Norsigian *et al*., (2010)](https://doi.org/10.1038/nprot.2009.203)        |     m9_spec.json                  |
| nutrient      |  This study. Peptone: [BD, 2015](https://legacy.bd.com/ds/technicalCenter/misc/lcn01558-bionutrients-manual.pdf), [Loginova *et al*., (1974)](https://doi.org/10.1007/BF00777001), [ThermoFisher, (2019)](https://assets.thermofisher.com/TFS-Assets/BPD/brochures/peptones-supplements-feeds-technical-reference-guide.pdf). Beef extract: [BD, 2015](https://legacy.bd.com/ds/technicalCenter/misc/lcn01558-bionutrients-manual.pdf), [ThermoFisher, (2019)](https://assets.thermofisher.com/TFS-Assets/BPD/brochures/peptones-supplements-feeds-technical-reference-guide.pdf)           |        nutrient_spec.json               |
|   pmm5_mendoza    |   [Mendoza *et al*., (2019)](https://doi.org/10.1186/s13059-019-1769-1)        |         pmm5_mendoza_spec.json              |
|  pmm7_mendoza     |  [Mendoza *et al*., (2019)](https://doi.org/10.1186/s13059-019-1769-1)         |    pmm7_mendoza_spec.json                   |
|  tsa     | This study. Tryptic soy: [BD, 2015](https://legacy.bd.com/ds/technicalCenter/misc/lcn01558-bionutrients-manual.pdf), [ThermoFisher, (2019)](https://assets.thermofisher.com/TFS-Assets/BPD/brochures/peptones-supplements-feeds-technical-reference-guide.pdf), [Hagely *et al*., (2013)](https://doi.org/10.1021/jf303985q), [Choct *et al*., (2010)](https://doi.org/10.5713/ajas.2010.90222)      |     tsa_spec.json                  |
| tsa_sheep_blood      |   This study. Tryptic soy: [BD, 2015](https://legacy.bd.com/ds/technicalCenter/misc/lcn01558-bionutrients-manual.pdf), [ThermoFisher, (2019)](https://assets.thermofisher.com/TFS-Assets/BPD/brochures/peptones-supplements-feeds-technical-reference-guide.pdf), [Hagely *et al*., (2013)](https://doi.org/10.1021/jf303985q), [Choct *et al*., (2010)](https://doi.org/10.5713/ajas.2010.90222)        |    tsa_sheep_blood_spec.json                   |


### Options

#### Required

`--assembly_fp` - Input assembly for which a metabolic model will be generated. This can be either an unannotated
**fasta** file or an annotated **genbank** file. Bactabolize will honour the genbank annotations. **IMPORTANT:** If you are using draft genome assemblies it is important to consider assembly quality and completeness. For *Klebsiella pneumoniae* we recommended minimum assembly quality: ≤200 assembly graph dead ends (calculate from .gfa or fastg). If only contigs are available, ≤130 contigs.

`--output_fp` - Output filename

`--ref_model_fp` - Reference model in .json format. (You can find the latest *K. pneumoniae* pan-model [here](https://github.com/kelwyres/KpSC-pan-metabolic-model).

`--biomass_reaction_id` - reaction ID of reference model's Biomass function. DEFAULT: BIOMASS_

The sequence data corresponding to the reference model can be provided as nucleotide (.ffn) AND protein multifasta (.faa) files (useful if the reference is a pan-model or multi-strain model for which the sequences do not exist in a single genbank file). Alternatively, the reference sequence data can be provided by a genbank file (useful if the reference is a single-strain model).

`--ref_genes_fp` - Reference genes in nucleotide sequence (fasta)

AND

`--ref_proteins_fp` - Reference protein in amino acid sequence (fasta).

OR

`ref_genbank_fp` - Reference genome (genbank).

#### Optional

`--media_type` - Growth media for initial model growth simulation test. One of: cdm_mendoza, bg11, lb, lb_carveme, m9, nutrient, pmm5_mendoza, pmm7_mendoza, tsa,
tsa_sheep_blood. DEFAULT: m9

`--atmosphere_type` - Choose atmosphere for model building. One of: aerobic, anaerobic. DEFAULT: aerobic

`--min_coverage` - Set minimum query coverage (%) for for ortholog identification via bi-directional best hit blast+. DEFAULT: 25

`--min_pident` - Set minimum identity (%) for for ortholog identification via bi-directional best hit blast+. DEFAULT: 80

`--min_ppos` - Set minimum protein similarity (positives, %) for for ortholog identification via bi-directional best hit blast+. DEFAULT: OFF. Can be used instead of `--min_pident` to allow for greater tolerance of similarly-functional but different residues.

`--memote_report_fp` - output file path for [MEMOTE model quality report](https://github.com/opencobra/memote). Note that this will significantly increase compute time e.g. +5 minutes PER assembly on a standard 1.60GHz laptop.

`--no_reannotation` - Will prevent the re-annotation of input genome assemblies if provided in genbank format. DEFAULT: Off. If set to `ON` [Prodigal](https://github.com/hyattpd/Prodigal) is used to identify coding regions.

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
    --output_fp input_assembly_qc_25_sim_85 \
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
    --output_fp input_assembly_qc_25_sim_85 \
    --min_coverage 25 \
    --min_pident 75 \
    --memote_report_fp input_report
```

### Draft model outputs

| Filename                      | Description                                      |
| ---------                     |---------                                         |
| `assembly_id.gbk`             | Prodigal annotation of input assembly            |
| `assembly_id_model.json`     | Metabolic model in .json format                  |
| `assembly_id_model.xml`      | Metabolic model in .xml ([SMBL Level 3 Version 1)](https://synonym.caltech.edu/documents/specifications/level-3/)) |
| `assembly_id_model.html`    | [MEMOTE](https://github.com/opencobra/memote) model report                             |

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

An `assembly_id_fba.tsv` tab-delimited file will be produced:

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
# K_variicola_variicola_342.json fails to produce biomass as it lacks DTDP-4-dehydrorhamnose 3,5-epimerase and
# DTDP-4-dehydrorhamnose reductase, which are required to create DTDP-L-rhamnose
# Here we assume DTDP-L-rhamnose is not essential for growth and patch the model accordingly
bactabolize patch_model \
    --draft_model_fp K_variicola_variicola_342.json \
    --ref_model_fp iYL1228_annotated.json \
    --biomass_reaction_id YOUR_BIOMASS_ID \
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

`--biomass_reaction_id` - ID of reference model's Biomass function. DEFAULT: BIOMASS_

##### Optional

`--media_type` - Choose growth media for model building. One of: bg11, lb_carveme, lb, m9, nutrient, tsa,
tsa_sheep_blood. DEFAULT: m9

`--atmosphere_type` - Choose atmosphere for model building. One of: aerobic, anaerobic. DEFAULT: aerobic

`--memote_report_fp` - MEMOTE model quality report output filepath. Note that this will add >5 minutes of compute time
PER assembly

## Metabolite IDs

To run individual FBA on extracellular metabolites, they must be annotated with the respective chemical formula in the
model. If you have a model that contains `metanetx` identifiers for metabolites (i.e. a BiGG model), you can add
metabolite formulas using the [`BiGG model compound
annotator`](https://github.com/scwatts/bigg_model_compound_annotator).

## Creation of custom medias
Creating custom medias depends on the Bactabolize module you are using. If you want to add a custom media for `draft_model` or `patch_model` or `sgk`, go [here](#draft-patch-or-sgk-custom-medias). For `fba`, go [here](#fba-custom-medias).

### Draft, patch or sgk custom medias

1. Copy m9 media file as template to your working directory

```
cp miniconda3/envs/bactabolize/lib/{your-python-version}/site-packages/bactabolize/data/media_definitions/m9_media.json your_working_dir/custom_media.json
```
2. Make edits to the `custom_media.json`. Be sure to edit the `name` as well as the `exchanges`.
3. Copy  `custom_media.json` into environment folder

```
cp your_working_dir/custom_media.json miniconda3/envs/bactabolize/lib/{your-python-version}/site-packages/bactabolize/data/media_definitions/
```
4. Run `bactabolize draft_model`, `patch_model` or `sgk` with `--media_type custom_media`


### FBA custom medias

1. Copy m9 spec file as template to your working directory

```
cp miniconda3/envs/bactabolize/lib/{your-python-version}/site-packages/bactabolize/data/fba_specs/m9_spec.json your_working_dir/custom_spec.json
```
2. Make edits to the custom_spec.json. Be sure to edit the `M9` dict header, along with `media_type`
3. Copy `custom_spec.json` into environment folder
```
cp your_working_dir/custom_spec.json miniconda3/envs/bactabolize/lib/{your-python-version}/site-packages/bactabolize/data/fba_specs/
```
4. Run `bactabolize fba` with `--fba_spec_fp data/fba_specs/custom_spec.json`


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
