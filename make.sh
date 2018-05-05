#!/bin/bash
#===============================================================================
#
#          FILE:  make.sh
# 
#         USAGE:  ./make.sh 
# 
#   DESCRIPTION:  create python extensions  
# 
#        AUTHOR:  Peter Robert Schmitt (iMac), pschmitt@upenn.edu
#       COMPANY:  University of Pennsylvania
#       VERSION:  1.0
#       CREATED:  06/11/2016 17:14:35 EDT
#      REVISION:  ---
#===============================================================================

echo "*** building MultiSURF "
python setup_multisurf.py build_ext -i
echo "*** building ReliefF "
python setup_relieff.py build_ext -i
echo "*** building SURF "
python setup_surf.py build_ext -i
