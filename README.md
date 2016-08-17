# Relief_Based_Algorithms
This suite includes stand-alone Python code to run any of
the included/available Relief-Based algorithms designed for attribute
filtering/ranking.  These algorithms are a quick way to identify
attributes in the dataset that may be most important to predicting some
phenotypic endpoint.  These scripts output an ordered set of attribute
names, along with respective scores (uniquely determined by the particular
algorithm selected).  Certain algorithms require key run parameters to
be specified.  

## MDR
This code is largely based on the Relief-Based algorithms
implemented in the Multifactor Dimensionality Reduction (MDR) software.
However these implementations have been expanded to accomodate continuous
attributes (and continuous attributes mixed with discrete attributes)
as well as a continuous endpoint.  This code also accomodates missing
data points.  Built into this code, is a strategy to automatically detect
from the loaded data, these relevant characteristics.

## Cython
Portions of this code has been optimized with Cython routines.

## How to run
1. Compile cython file(s)  ./make.sh
2. Run rebate.py -h for all the available options
3. No Windows support as yet.

## To Do
1. Make this an installable package
2. Convert to Classes
3. Create GUI Interface
