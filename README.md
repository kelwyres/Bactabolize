<img src="https://user-images.githubusercontent.com/19924405/193505313-edd9453a-e4eb-4730-81b1-a2bd9e652721.png" width="50%">

A high-throughput genome-scale metabolic reconstruction and growth simulation pipeline.

Bactabolize is designed for rapid generation of strain-specific metabolic reconstructions from bacterial genome data using the approach described in [Norsigian et al. Nature Protocols 2020](https://www.nature.com/articles/s41596-019-0254-3). It leverages the [COBRApy toolkit](https://opencobra.github.io/cobrapy/) and takes an input genome assembly (annotated or unannotated) to construct a strain-specific draft model by comparison to a reference (ideally a multi-strain or 'pan'-reference model). It also allows high-throughput growth phenotype simulation via Flux Balance Analysis (FBA) e.g. to predict substrate usage profiles and Single Gene Knockout analysis (SGK) to predict the impacts of single gene knockout mutations. These can be performed under a variety of growth conditions and mediums.

Bactabolize is freely available under a [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html). Please cite the following papers if you make use of Bactabolize:

* Vezina B. / Watts S.C. et al. Bactabolize. _In prep_
* Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC
  Syst Biol 7, 74 (2013). <https://doi.org/10.1186/1752-0509-7-74>

**Visit the wiki to find out more!**
