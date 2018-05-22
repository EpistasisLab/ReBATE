Master status: [![Master Build Status](https://travis-ci.org/EpistasisLab/ReBATE.svg?branch=master)](https://travis-ci.org/EpistasisLab/ReBATE)
[![Master Coverage Status](https://coveralls.io/repos/github/EpistasisLab/ReBATE/badge.svg?branch=master)](https://coveralls.io/github/EpistasisLab/ReBATE?branch=master)

Development status: [![Development Build Status](https://travis-ci.org/EpistasisLab/ReBATE.svg?branch=development)](https://travis-ci.org/EpistasisLab/ReBATE)
[![Development Coverage Status](https://coveralls.io/repos/github/EpistasisLab/ReBATE/badge.svg?branch=development&service=github)](https://coveralls.io/github/EpistasisLab/ReBATE?branch=development)

Package information: ![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)
![Python 3.5](https://img.shields.io/badge/python-3.6-blue.svg)
![License](https://img.shields.io/badge/license-MIT%20License-blue.svg)

# ReBATE (Relief-based Algorithm Training Environment) 
This package includes stand-alone Python code to run any of the included/available [Relief-Based algorithms (RBAs)](https://en.wikipedia.org/wiki/Relief_(feature_selection)) 
designed for feature weighting/selection as part of a machine learning pipeline (supervised learning). Presently this includes the following core RBAs: ReliefF, SURF, SURF\*, MultiSURF\* and MultiSURF. Additionally, an implementation of the iterative TuRF mechanism is included. **It is still under active development** and we encourage you to check back on this repository regularly for updates.

These algorithms offer a computationally efficient way to perform feature selection that is sensitive to feature interactions as well as simple univariate associations, unlike most currently available filter-based feature selection methods. The main benefit of Relief algorithms is that they identify feature interactions without having to exhaustively check every pairwise interaction, thus taking significantly less time than exhaustive pairwise search.

Each core algorithm outputs an ordered set of feature names along with respective feature scores (i.e. weights). Certain algorithms require user specified run parameters (e.g. ReliefF requires the user to specify some 'k' number of nearest neighbors). 

Relief algorithms are commonly applied to genetic analyses, where epistasis (i.e., feature interactions) is common. However, the algorithms implemented in this package can be applied to almost any supervised classification data set and supports:

* Feature sets that are discrete/categorical, continuous-valued or a mix of both

* Data with missing values

* Binary endpoints (i.e., classification)

* Multi-class endpoints (i.e., classification)

* Continuous endpoints (i.e., regression)

Built into this code, is a strategy to 'automatically' detect from the loaded data, these relevant characteristics.

Of our two initial ReBATE software releases, this stand-alone version primarily focuses on improving run-time with the use of Cython. 
This code is most appropriate for more experienced users or those primarily interested in reducing analysis run time. 

We recommend that scikit-learn users, Windows operating system users, beginners, or those looking for the most recent ReBATE developments to instead use our alternate [scikit-rebate](https://github.com/EpistasisLab/scikit-rebate) implementation. ReBATE can be run on Windows with some additional installation steps and possible troubleshooting outlined below. 

## License
Please see the [repository license](https://github.com/EpistasisLab/ReBATE/blob/master/LICENSE) for the licensing and usage information for ReBATE.
Generally, we have licensed ReBATE to make it as widely usable as possible.

## Cython (Important Notice)
NOTICE: As is, this code will not run on your local platform! Portions of this code have been optimized with Cython routines for code speedup. As a result, before being able to use 
ReBATE on a given operating system (i.e. Linux, Mac, or Windows), critical binary files must be compiled as a one time step (or any time the underlying source code is modified, or
any time an updated version of ReBATE is downloaded to your system.  Compiling the necessary binary files is very easy to do on Mac or Linux systems (because they include a C compiler). However Windows users will unfortunately have to go through a few extra hurdles in order to complete this one time step. If you wish to avoid this hassle, please see our alternate [scikit-rebate](https://github.com/EpistasisLab/scikit-rebate) implementation.

## Installation
For detailed information on installing ReBATE, including necessary prerequisites, special instructions for Windows users, and instructions for compiling cython, please refer to our [installation documentation](https://EpistasisLab.github.io/ReBATE/installing/). 

## Running ReBATE
From the '/rebate/' directory, run the following to view all available options: 
```
./rebate.py -h
```

For detailed information and examples of how to run the different Relief algorithms available in this package, please refer to our [usage documentation](https://EpistasisLab.github.io/ReBATE/using/).

## Contributing to ReBATE
We welcome you to [check the existing issues](https://github.com/EpistasisLab/ReBATE/issues/) for bugs or enhancements to work on. If you have an idea for an extension to ReBATE, please [file a new issue](https://github.com/EpistasisLab/ReBATE/issues/new) so we can discuss it.

If you wish to contribute to ReBATE we strongly recommend following the steps detailed in [contributing documentation](https://EpistasisLab.github.io/ReBATE/contributing/).

## Citing ReBATE

If you use ReBATE or the MultiSURF algorithm in a scientific publication, please consider citing the following paper:

Ryan J. Urbanowicz, Randal S. Olson, Peter Schmitt, Melissa Meeker, Jason H. Moore (2017). [Benchmarking Relief-Based Feature Selection Methods](https://arxiv.org/abs/1711.08477). *arXiv preprint*, under review.

BibTeX entry:

```bibtex
@misc{Urbanowicz2017Benchmarking,
    author = {Urbanowicz, Ryan J. and Olson, Randal S. and Schmitt, Peter and Meeker, Melissa and Moore, Jason H.},
    title = {Benchmarking Relief-Based Feature Selection Methods},
    year = {2017},
    howpublished = {arXiv e-print. https://arxiv.org/abs/1711.08477},
}
```
If you wish to directly cite the original paper for one of the other algorithms implemented in ReBATE please refer to our [citing documentation](https://epistasislab.github.io/ReBATE/citing/).


## History
This code is largely based on Python implementations of ReliefF, SURF, SURF*, MultiSURF*, and TuRF within the [ExSTraCS](https://github.com/ryanurbs/ExSTraCS_2.0) algorithm software. 
That Python code was in turn based on Java implementations of these algorithms within the [Multifactor Dimensionality Reduction (MDR)](https://sourceforge.net/projects/mdr/) software.
In contrast with the MDR implementations, both the ExSTraCS and scikit-rebate, and present ReBATE versions of this code have been expanded to accommodate the following data considerations: Continuous features, a mix of discrete and continuous features, a continuous endpoint/outcome, and missing data values.

## Possible future updates
1. Make this an installable package
2. Convert to Classes
3. Create GUI Interface
