<img src="https://user-images.githubusercontent.com/19924405/193505313-edd9453a-e4eb-4730-81b1-a2bd9e652721.png" width="50%">

A high-throughput genome-scale metabolic reconstruction and growth simulation pipeline.

## How to run
**Install and quick start [here](https://github.com/kelwyres/Bactabolize/wiki/1.-Quick-start)**

**Visit the [wiki](https://github.com/kelwyres/Bactabolize/wiki) to find out more!**

## Description
Bactabolize is designed for rapid generation of strain-specific metabolic reconstructions from bacterial genome data
using the approach described in [Norsigian et al. Nature Protocols
2020](https://www.nature.com/articles/s41596-019-0254-3). It leverages the [COBRApy
toolkit](https://opencobra.github.io/cobrapy/) and takes an input genome assembly (annotated or unannotated) to
construct a strain-specific draft model by comparison to a reference (ideally a multi-strain or 'pan'-reference model).
It also allows high-throughput growth phenotype simulation via Flux Balance Analysis (FBA) e.g. to predict substrate
usage profiles and Single Gene Knockout analysis (SGK) to predict the impacts of single gene knockout mutations. These
can be performed under a variety of growth conditions and mediums.

## Compatible pan-reference metabolic models

|Species                        |Database                                                      |Reference    |
|-----------                    |-----------                                                   |-----------  |
|_Klebsiella pneumoniae_ Species Complex|[Model & Associated Sequences](https://github.com/kelwyres/KpSC-pan-metabolic-model)|[Cooper 2024, _MGen_](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001206)|
|_Bacillus subtilis_            |[Model](https://github.com/SBGlab/Bacillus_Subtilis_multistrain_GEM/tree/main/iBB1018) & [Gene Annotations](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009045.1/)  |[Blázquez 2023, _Int. J. Mol. Sci._](https://www.mdpi.com/1422-0067/24/8/7091)|

## Licence and citation

Bactabolize is freely available under a [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).
Please cite the following papers if you make use of Bactabolize:

* Vezina B. / Watts S.C. et al. 'Bactabolize: A tool for high-throughput generation of bacterial strain-specific metabolic models'. eLife (2023). 
  [https://doi.org/10.7554/eLife.87406.3](https://doi.org/10.7554/eLife.87406.3)
* Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. 'COBRApy: COnstraints-Based Reconstruction and Analysis for Python'. BMC
  Syst Biol 7, 74 (2013). <https://doi.org/10.1186/1752-0509-7-74>

