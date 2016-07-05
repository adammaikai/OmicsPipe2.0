#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages
#from distutils.core import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
  name = 'omics_pipe',
  version = '1.1.3', #major.minor[.patch[.sub]]
  author = 'Kathleen Fisch, Ph.D.',
  author_email = 'kfisch@scripps.edu',
  description = 'Extensible computational pipeline for running next generation sequencing analyses',
  license = 'MIT',
  keywords = "biology omics transcriptomics pipeline RNA-seq ChIP-seq WGS WES genome exome HPC",
  url = 'https://bitbucket.org/sulab/omics_pipe',
  download_url = "http://pypi.python.org/pypi/omics_pipe",
  long_description="Omics pipe is an open-source, modular computational platform that automates best practice multi-omics data analysis pipelines published in Nature Protocols and other commonly used pipelines, such as GATK.  It currently automates and provides summary reports for two RNAseq and miRNAseq pipelines, variant calling from whole exome sequencing (WES), variant calling and copy  number variation analysis from whole genome sequencing (WGS), two ChIPseq pipelines and a custom RNAseq pipeline for personalized genomic medicine reporting.  It also provides automated support for interacting with the The Cancer Genome Atlas (TCGA) datasets, including automatic download and processing of the samples in this database.",
  packages = ['omics_pipe'],
  #package_data = {'omics_pipe' : ['scripts/reporting/data/brca_mol_class/*','scripts/reporting/data/DoG/*', 'scripts/reporting/data/geneLists/*','scripts/reporting/ref/*', 'scripts/reporting/src/*']},
  include_package_data = True, #include everything in source control
  entry_points= {
        'console_scripts': [
            'omics_pipe = omics_pipe.main:main'
                            ]
                 },
  install_requires=[
    "ruffus>=2.0.0",
    "drmaa",
    "sumatra>=0.5.2",
    "PyYAML>=3.10",
    "argparse",
    "ez_setup",
    "setuptools",
    "mercurial",
    "hgapi",
    "raven",
    "importlib"
    
  ],
  classifiers = [
    'Programming Language :: Python',
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: MIT License',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    
  ],
)

