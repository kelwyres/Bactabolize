# Install python package
${PYTHON} -m pip install -vv .

# Install MEMOTE with pip; no Conda package available.
# NOTE(SW): several pip environment variables are set by conda-build, revert effective use of
# --no-index to allow discovery of the MEMOTE package on pypi. For more details see:
#  * https://pip.pypa.io/en/stable/topics/configuration/#environment-variables
#  * https://github.com/conda/conda-build/blob/e4d9b3/conda_build/build.py#L2593
PIP_NO_INDEX=False ${PYTHON} -m pip install 'memote ==0.13.0'
