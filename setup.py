#!/usr/bin/env python
# -*- coding: utf-8 -*-


# Always prefer setuptools over distutils
from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path
import subprocess
import sys
import time
import shutil
import os
		


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.txt'), encoding='utf-8') as f:
	long_description = f.read()

setup(
    name='paneltimegui',
    version='1.1.16',
    description='An efficient integrated panel and GARCH estimator',
    long_description=long_description,
    url='https://github.com/espensirnes/paneltimegui',
    author='Espen Sirnes',
    author_email='espen.sirnes@uit.no',
    license='GPL-3.0',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Researchers',
        'Topic :: Statistical Software :: time series estimation',
        'License :: OSI Approved :: GPL-3.0 License',
        'Programming Language :: Python :: 3.5',
        ],

    keywords='econometrics',
    packages=find_packages(exclude=[]),
    install_requires=['numpy ','tk','matplotlib','pymysql'],
	extras_require={'linux':'gcc'},	

    package_data={
        '': [],
        },
    include_package_data=True,

    entry_points={
        'console_scripts': [
            'paneltimegui=paneltimegui:main',
            ],
        },
)
