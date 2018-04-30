#!/bin/bash
#===============================================================================
#
#          FILE:  clean.sh
# 
#         USAGE:  ./clean.sh 
# 
#   DESCRIPTION:  Clean up for recompiling
# 
#        AUTHOR:  Peter Robert Schmitt (@HOME), pschmitt@upenn.edu
#       COMPANY:  University of Pennsylvania
#       VERSION:  1.0
#       CREATED:  05/23/2016 07:24:26 EDT
#===============================================================================

rm -rvf build *.c *.so *.pyc __pycache__ data/*scores*
