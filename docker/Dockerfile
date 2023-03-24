FROM mambaorg/micromamba:0.24.0

RUN \
  micromamba install -y -n base -c scwatts -c bioconda -c conda-forge -c defaults \
    'bactabolize==1.0.1' && \
  micromamba clean --all --yes

ENV PATH="/opt/conda/bin:/opt/conda/condabin:${PATH}"
