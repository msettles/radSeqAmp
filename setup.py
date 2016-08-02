#!/usr/bin/env python

'''
radSeqAmp
Analysis of radSeq data, Restriction Associated Digest Sequecnign, initial
version is to process the samples (trim) barcodes, detect restriction site

Later version will exand on analysis
'''

import sys

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension


editdist = Extension('editdist', sources=['lib/editdist.c'])
trim = Extension('trim', sources=['lib/trim.c'])

try:
    version_num = open("VERSION", "r+").readline().strip()
except:
    sys.stderr.write("Error retrieving version_number")


config = \
    {
        'description': 'Processing of Illumina amplicon projects - radSeq version',
        'author': 'Matt Settles',
        'url': 'https://github.com/msettles/radSeqAmp',
        'download_url': 'https://github.com/msettles/radSeqAmp',
        'author_email': 'settles@ucdavis.edu',
        'version': version_num,
        'install_requires': [],
        'packages': ['radSeqAmp'],
        'scripts': ['bin/radSeqAmp'],
        'name': 'radSeqAmp',
        "ext_package": 'radSeqAmp',
        'ext_modules': [editdist, trim]
    }

setup(**config)
