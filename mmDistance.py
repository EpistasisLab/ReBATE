# Wed Aug 24 14:51:56 EDT 2016
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
