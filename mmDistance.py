# Wed Aug 24 14:51:56 EDT 2016
"""
Copyright (c) 2016 Peter R. Schmitt and Ryan J. Urbanowicz

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
"""
import sys
import numpy as np
from numpy import isnan, where, append, unique, delete, empty
###############################################################################
def getDistances(xc, xd, var, cdiffs):
    """ This creates a distance array for mixed types of data with or 
        without missing data """

    distArray = []
    datalen = var['datalen']
    missing = int(var['mdcnt'])
    
    # get indices of missing data per record
    if(missing > 0):
        cindices = list()
        dindices = list()
        for i in range(datalen):
            cindices.append(where(isnan(xc[i]))[0])
            dindices.append(where(isnan(xd[i]))[0])

    
    for index in range(datalen):

        if(missing > 0):
            row = getrow_missing(xc, xd, cdiffs, index, cindices, dindices)
        else:
            row = getrow_mixed(xc, xd, cdiffs, index)

        row = list(row)
        distArray.append(row)
        
    return distArray
###############################################################################
def getrow_missing(xc, xd, cdiffs, index, cindices, dindices):

    row = empty(0,dtype=np.double)
    cinst1 = xc[index]
    dinst1 = xd[index]
    can = cindices[index]
    dan = dindices[index]
    for j in range(index):
        dist = 0
        dinst2 = xd[j]
        cinst2 = xc[j]

        # continuous
        cbn = cindices[j]
        idx = unique(append(can,cbn))   # create unique list
        c1 = delete(cinst1,idx)       # remove elements by idx
        c2 = delete(cinst2,idx)
        cdf = delete(cdiffs,idx)

        # discrete
        dbn = dindices[j]
        idx = unique(append(dan,dbn))
        d1 = delete(dinst1,idx)
        d2 = delete(dinst2,idx)
            
        # discrete first
        dist += len(d1[d1 != d2])

        # now continuous
        dist += np.sum(np.absolute(np.subtract(c1,c2)) / cdf)

        row = append(row,dist)

    return row
##############################################################################
# if data is clean this will not be used
def getrow_mixed(xc, xd, cdiffs, index):

    row = empty(0,dtype=np.double)
    d1 = xd[index]
    c1 = xc[index]
    for j in range(index):
        dist = 0
        d2 = xd[j]
        c2 = xc[j]

        # discrete first
        dist += len(d1[d1 != d2])

        # now continuous
        dist += np.sum(np.absolute(np.subtract(c1,c2)) / cdiffs)

        row = append(row,dist)

    return row
##############################################################################
