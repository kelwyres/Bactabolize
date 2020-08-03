# metabolic
A high-throughput metabolic model construction pipeline.


## Table of contents
* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [License](#license)


## Quickstart
```bash
# Clone the repository
git clone https://github.com/scwatts/metabolic_pipeline.git && cd metabolic_pipeline/

# Install dependencies and activate conda environment
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file misc/conda_packages.txt
conda activate $(pwd -P)/conda_env

# Run
./metabolic-runner.py
```

## Requirements
### Reference model
To run individual FBA on extracellular metabolites, they must be annotated with the respective chemical formula in the model.
If you have a model that contains `metanetx` identifiers for metabolites (i.e. a BiGG model), you can add metabolite formulas
using the [`BiGG model compound annotator`](https://github.com/scwatts/bigg_model_compound_annotator).


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
