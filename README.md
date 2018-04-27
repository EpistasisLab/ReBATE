Package information: ![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)
![Python 3.5](https://img.shields.io/badge/python-3.6-blue.svg)
![License](https://img.shields.io/badge/license-MIT%20License-blue.svg)

# ReBATE (Relief-based Algorithm Training Environment) 
This package includes stand-alone Python code to run any of the included/available [Relief-Based algorithms (RBAs)](https://en.wikipedia.org/wiki/Relief_(feature_selection)) 
designed for feature weighting/selection as part of a machine learning pipeline (supervised learning). Presently this includes the following core RBAs: ReliefF, SURF, SURF*, and MultiSURF*. Additionally, an implementation of the iterative TuRF mechanism is included. **It is still under active development** and we encourage you to check back on this repository regularly for updates.

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
### Prerequisites
All of the necessary Python packages can be installed via the [Anaconda Python distribution](https://www.continuum.io/downloads), which we strongly recommend that you use. We also strongly recommend that you use Python 3 over Python 2 if you're given the choice.

ReBATE requires that the following external Python packages be installed (all included in Anaconda):

argparse, time, sys, os, IO, Common, numpy, math, pandas, scipy.spatial.distance, operator, csv, distutils.core, distutils.extension, Cython.Distutils, datetime

NumPy and SciPy can be installed in Anaconda via the command:

```
conda install numpy scipy
```
#### Additional Steps For Windows Users
Here we discuss the additional prerequisites for running ReBATE on a Windows operating system. First you will need a command line terminal for Windows (e.g. [Cygwin](https://www.cygwin.com/install.html), or [GitBash](https://gitforwindows.org/)) with Anaconda properly installed. From this terminal you can compile the Cython code, and run ReBATE. Second you will need a C compiler for compiling the Cython code once before being able to run ReBATE on your Windows machine. 

There are a number of possible ways to get the C compiler working on your Windows terminal, all of which will depend on the version of windows, the version of python/Anaconda, and whether it is 32-bit or 64-bit. Be aware that there may be some troubleshooting in getting the C compiler operational in your command line terminal. Below we outline steps that worked for us in the spring of 2018 using Windows 10, and Python 3.5.2 with Anaconda 4.0 (64-bit), using GitBash as our terminal.

1.) Ensure that setuptools is updated. Run the following in your terminal: pip install -upgrade setuptools

2.) [Download/install Visual Studio Community 2017](https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2017). This includes the necessary C compiler. It is not necessary to install all of the visual studio components, as this can be slow and take up a good deal of space. However make sure that you download and install the 'Desktop Development with C++' and the 'Python Development' workloads as detailed at the following [link](https://docs.microsoft.com/en-us/visualstudio/python/working-with-c-cpp-python-in-visual-studio).  Within the Python Development workload, also select the box on the right for Python native development tools. 

3.) Once Visual Studio has been successfully installed, there was a remaining bug that required some troubleshooting for the C compiler to be found by the terminal. The following fix worked for us: 

(A) [Add this to your PATH environment variables](https://www.java.com/en/download/help/path.xml): C:\Program Files (x86)\Windows Kits\10\bin\x64   

(B) Copy these two files (rc.exe and rcdll.dll) from (C:\Program Files (x86)\Windows Kits\8.1\bin\x86) to (C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin)

4.) At this point you should be able to compile Cython as described below.

### Compile Cython
Once these prerequisites are installed, it will be necessary to compile the Cython code on the respective operating system within which ReBATE will be run. It is only necessary to do this once, not every time ReBATE is run. This happens in two stages (1) a .pyx file is compiled by cython to a .c file, then (2) the .c file is compiled by a C compiler to a .so file (or a .pyd file for Windows). 

Simply run the following file included with ReBATE to produce the .c and (.so or .pyd) files:  ./make.sh

If there is need to recompile the Cython files, first remove the previous .c and (.so or .pyd) files by running: ./clean.sh

## Contributing to ReBATE

We welcome you to [check the existing issues](https://github.com/EpistasisLab/ReBATE/issues/) for bugs or enhancements to work on. If you have an idea for an extension to ReBATE, please [file a new issue](https://github.com/EpistasisLab/ReBATE/issues/new) so we can discuss it.


## Citing ReBATE

If you use ReBATE in a scientific publication, please consider citing the following paper:

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

## History
This code is largely based on Python implementations of ReliefF, SURF, SURF*, MultiSURF*, and TuRF within the [ExSTraCS](https://github.com/ryanurbs/ExSTraCS_2.0) algorithm software. 
That Python code was in turn based on Java implementations of these algorithms within the [Multifactor Dimensionality Reduction (MDR)](https://sourceforge.net/projects/mdr/) software.
In contrast with the MDR implementations, both the ExSTraCS and scikit-rebate, and present ReBATE versions of this code have been expanded to accommodate the following data considerations: Continuous features, a mix of discrete and continuous features, a continuous endpoint/outcome, and missing data values.

## Possible future updates
1. Make this an installable package
2. Convert to Classes
3. Create GUI Interface
