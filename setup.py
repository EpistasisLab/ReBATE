#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

def calculate_version():
    initpy = open('_version.py').read().split('\n')
    version = list(filter(lambda x: '__version__' in x, initpy))[0].split('\'')[1]
    return version

package_version = calculate_version()

setup(
    name='ReBATE',
    version=package_version,
    author='Pete Schmitt, and Ryan J. Urbanowicz',
    author_email='ryanurb@upenn.edu',
    packages=find_packages(),
    url='https://github.com/EpistasisLab/ReBATE',
    license='License :: OSI Approved :: MIT License',
    description=('Relief-based feature selection algorithms'),
    long_description='''
A Cython optimized Python implementation of ReBATE, a suite of Relief-based feature selection algorithms.

Contact
=============
If you have any questions or comments about ReBATE, please feel free to contact us via e-mail: ryanurb@upenn.edu

This project is hosted at https://github.com/EpistasisLab/ReBATE
''',
    zip_safe=True,
    install_requires=['numpy', 'scipy', 'scikit-learn'],
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Utilities'
    ],
    keywords=['data mining', 'feature selection', 'feature importance', 'machine learning', 'data analysis', 'data engineering', 'data science'],
    include_package_data=True,
)
