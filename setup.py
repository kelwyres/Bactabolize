#!/usr/bin/env python3
import setuptools
import sys


import metabolic


setuptools.setup(
    name='metabolic',
    version=metabolic.__version__,
    description='metabolic python package',
    author='Stephen Watts',
    license='GPLv3',
    test_suite='tests',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': ['metabolic=metabolic.__main__:entry'],
    }
)
