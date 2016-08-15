# No longer used
from __future__ import division
import numpy as np
from math import isnan
#######################################
# use fabs function from c math library
#######################################
cdef extern from "math.h":
    double fabs(double)

cdef double get_abs(double value):
    return fabs(value)
#######################################

###############################################################################
cdef double compare_mixed(double[::1] inst1, double[::1] inst2, 
                             int[::1] dtypes, double[::1] diffs):
    cdef:
        double d = 0.0
        int ilen = len(inst1)
        int i

    for i in range(ilen):
        if(dtypes[i] == 0):  # dtype is discrete
            if(int(inst1[i]) != int(inst2[i])): d += 1.0
        else:                # dtype is continuous
            d += get_abs(inst1[i] - inst2[i]) / diffs[i]

    d /= ilen
    return d
###############################################################################
cdef double compare_missing(double[::1] inst1, double[::1] inst2, 
                            int[::1] dtypes, double[::1] diffs):
    cdef:
        double d = 0.0
        int ilen = len(inst1)
        int nancnt = 0
        int i

    for i in range(ilen):
        if(dtypes[i] == 0):  # dtype is discrete
            if(isnan(inst1[i]) or isnan(inst2[i])):
                nancnt += 1
            elif(int(inst1[i]) != int(inst2[i])): 
                d += 1.0
        else:                # dtype is continuous
            if(isnan(inst1[i]) or isnan(inst2[i])):
                nancnt += 1
            else:
                d += get_abs(inst1[i] - inst2[i]) / diffs[i]

    d = d / (ilen - nancnt)
    return d 

###############################################################################
# This creates a row distance array for all types of data
def get_row_dist(double[:,::1] mdata, int[::1] dtypes, double[::1] diffs, 
                 int index, int missing):

    dtypes = np.array(dtypes, dtype=np.int32)
    dlength = len(mdata)
    cdef int j
    cdef double dist

    row = []
    inst1 = mdata[index]
    if(missing > 0):
        for j in range(index):
            inst2 = mdata[j]
            dist = compare_missing(inst1, inst2, dtypes, diffs)
            row.append(dist)
    else:
        for j in range(index):
            inst2 = mdata[j]
            dist = compare_mixed(inst1, inst2, dtypes, diffs)
            row.append(dist)

    return row
###############################################################################
