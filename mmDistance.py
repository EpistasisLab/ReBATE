import sys
import numpy as np
from numpy import isnan, where, append, unique, delete, empty
from IO import printf
###############################################################################
def getDistances(xc, xd, var, cdiffs):
    """ This creates a distance array for mixed types of data with or 
        without missing data """

    distArray = []

    missing = int(var['mdcnt'])
    
    for index in range(var['datalen']):

        row = getrow(xd, xc, cdiffs, index, missing)
        row = list(row)
        distArray.append(row)
        
    return distArray
###############################################################################
def getrow(xd, xc, cdiffs, index, missing):
    
    row = empty(0,dtype=np.double)
    dinst1 = xd[index]
    cinst1 = xc[index]
    for j in range(index):
        dist = 0
        dinst2 = xd[j]
        cinst2 = xc[j]

        if(missing > 0):
            # continuous
            an = where(isnan(cinst1))[0]  # find idx of missing data
            bn = where(isnan(cinst2))[0]
            idx = unique(append(an,bn))   # create unique list
            c1 = delete(cinst1,idx)       # remove elements by idx
            c2 = delete(cinst2,idx)
            cdf = delete(cdiffs,idx)

            # discrete
            an = where(isnan(dinst1))[0]
            bn = where(isnan(dinst2))[0]
            idx = unique(append(an,bn))
            d1 = delete(dinst1,idx)
            d2 = delete(dinst2,idx)
        else:
            d1 = dinst1
            d2 = dinst2
            c1 = cinst1
            c2 = cinst2
            
        # discrete first
        dist += len(d1[d1 != d2])

        # now continuous
        dist += np.sum(np.absolute(np.subtract(c1,c2)) / cdf)

        row = append(row,dist)

    return row
##############################################################################
