#!/usr/bin/env python3
import setuptools


import bactabolize


setuptools.setup(
    name='Bactabolize',
    version=bactabolize.__version__,
    description='Bactabolize python package',
    author='Stephen Watts',
    license='GPLv3',
    url='https://github.com/scwatts/bactabolize',
    test_suite='tests',
    packages=setuptools.find_packages(),
    package_data={'bactabolize': ['data/*.json']},
    entry_points={
        'console_scripts': ['bactabolize=bactabolize.__main__:entry'],
    },
)
