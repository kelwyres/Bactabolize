#!/usr/bin/env python3
import setuptools


import metabolic


setuptools.setup(
    name='metabolic',
    version=metabolic.__version__,
    description='metabolic python package',
    author='Stephen Watts',
    license='GPLv3',
    url='https://github.com/scwatts/metabolic_pipeline',
    test_suite='tests',
    packages=setuptools.find_packages(),
    package_data={'metabolic': ['data/*.json']},
    entry_points={
        'console_scripts': ['metabolic=metabolic.__main__:entry'],
    },
)
