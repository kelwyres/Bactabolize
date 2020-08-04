# metabolic
A high-throughput metabolic model construction pipeline.


## Table of contents
* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [Requirements](#requirements)
* [Outputs](#outputs)
* [License](#license)


## Quickstart
```bash
# Clone the repository
git clone https://github.com/scwatts/metabolic_pipeline.git && cd metabolic_pipeline/

# Install dependencies and activate conda environment
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file misc/conda_packages.txt
conda activate $(pwd -P)/conda_env

# Run an example
./metabolic-runner.py --assembly_fp data/K_variicola_tropicalensis_CDC4241-71.gbk \
    --ref_genbank_fp data/K_pneumoniae_MGH78578.gbk \
    --ref_model_fp data/model_annotated.json \
    --fba_spec_fp data/m9_fba_spec.json \
    --prodigal_model_fp prodigal_models/output/klebsiella_pneumoniae_GCA_006364295.1_ASM636429v1.fasta_model.bin \
    --output_dir output/
```

## Requirements
### Reference model
To run individual FBA on extracellular metabolites, they must be annotated with the respective chemical formula in the model.
If you have a model that contains `metanetx` identifiers for metabolites (i.e. a BiGG model), you can add metabolite formulas
using the [`BiGG model compound annotator`](https://github.com/scwatts/bigg_model_compound_annotator).


## Outputs
| Filename                      | Description                           |
| ---------                     |---------                              |
| `assembly_id`\_stats.tsv      | Input assembly statistics             |
| `assembly_id`.gbk             | Prodigal annotation of input assembly |
| `assembly_id`\_model.json     | Metabolic model                       |
| `assembly_id`\_fba.tsv        | Results of FBA                        |


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
