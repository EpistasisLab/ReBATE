# ReBATE 0.2

* Added documentation with mkdocs.

* Added unit testing. Adopted Travis CI

* Added 'MultiSURF' algorithm, previously only available in scikit-rebate.

* Updated ReBATE to include other updates made to scikit-rebate. 

* Fixed score normalizations so they fall between -1 and 1 for all algorithms Now matches scikit-rebate.

* Consolidated MultiSURF* so that one script is used for both multiclass, and other types of endpoints. 

* Added an automatic (standard deviation based) ramp function method that is utilized by all algorithms on data with a mix of discrete and continuous features. Taken from scikit-rebate.

* Included steps to support operation of ReBATE in Windows. 

* Beyond what was previously used for testing scikit-rebate, we added the 6-bit multiplexer as a test problem as well as a simple 3-class (multiclass) GWAS-style simulated dataset (with 100% heritability). 

# ReBATE 0.1

* Initial release of Relief algorithms, including ReliefF, SURF, SURF*, MultSURF*, and TuRF.
