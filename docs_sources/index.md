[ReBATE](https://github.com/EpistasisLab/ReBATE) includes stand-alone Python code to run any of the included/available [Relief-Based algorithms (RBAs)](https://en.wikipedia.org/wiki/Relief_(feature_selection)) 
designed for feature weighting/selection as part of a machine learning pipeline (supervised learning). Presently this includes the following core RBAs: ReliefF, SURF, SURF\*, and MultiSURF\*. Additionally, an implementation of the iterative TuRF mechanism is included. As of 5/7/18, **it is still under active development** and we encourage you to check back on this repository regularly for updates.

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