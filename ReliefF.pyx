# Thu Jul 14 12:21:00 EDT 2016
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
import time as tm
import numpy as np
import sys
import math
###############################################################################
def runReliefF(header,double[:,::1] x,double[::1] y,attr,var,distArray,options):
    """This first finds the number of nearest neighbors specified by the value 
       of numNeighbors in the variables dictionary"""
    V = options['verbose']
    start = tm.time()
    #-------------------------------------------------------------------------
    def buildIndexLists():
        """ This creates lists of indexes of observations that share the same
            value in the phenotype"""
        cdef int i, index = 0
        indicies = dict()

        # initialize dictionary to hold as many lists
        # as there are unique values in the phenotype
        for i in var['phenoTypeList']:
            indicies[i] = []

        for i in y:
            indicies[i].append(index)
            index += 1

        return indicies
    #-------------------------------------------------------------------------
    def getdistance(i):
        if(i == inst):
            return sys.maxsize
        elif(i < inst):
            return distArray[inst][i]
        else:
            return distArray[i][inst]
    #-------------------------------------------------------------------------
    def getsortTuple(x):
        return (getdistance(x),x)
    #-------------------------------------------------------------------------
    # Find number of classes in the dataset and store them into the map
    def getMultiClassMap():
        mcmap = dict()

        for i in range(maxInst):
            if(y[i] not in mcmap):
                mcmap[y[i]] = 0
            else:
                mcmap[y[i]] += 1

        for each in var['phenoTypeList']:
            mcmap[each] = mcmap[each]/float(maxInst)
            
        return mcmap
    #-------------------------------------------------------------------------
    cdef:
        int neighbors = var['numNeighbors']
        int numattr = var['NumAttributes']
        int datalen = var['datalen']
        int maxInst = datalen
        
    Scores = [0] * numattr

    indicies = buildIndexLists()
    
    if(var['classType'] == 'multiclass'):
        mcmap = getMultiClassMap()
    else:
        mcmap = 0

    for inst in range(datalen):
        NN = []
        for sample in indicies.values():
            n = sorted(sample,key=getsortTuple)
            NN.extend(n[:neighbors])
        
        for ai in range(numattr):
            #nn = np.array(NN)
            nn = np.ascontiguousarray(NN, dtype=np.int32)
            Scores[ai] += getReliefFScores(header,attr,var,x,y,nn,inst,ai,mcmap)

    #averaging the scores
    divisor = maxInst * neighbors
    for ai in range(numattr):
        Scores[ai] = Scores[ai]/float(divisor)

    if(V):
        print("getNeighbors + getScores elapsed time = " \
               + str(tm.time() - start))
        sys.stdout.flush()
    return Scores
###############################################################################
#  Method evaluating the score of an attribute
#  called from runReliefF()
cdef double getReliefFScores(header, attr, var, double[:,::1] x, double[::1] y,
                             int[::1] NN, int inst, int ai, mcmap):

    cdef:
        double hit_count, miss_count, hit_diff, miss_diff, diff
        double missSum, missAvg, mmdiff, calc = 0
        double hit_prop, miss_prop 
        double same_class_bound = var['phenSD']
        double inst_item, NN_item
        int i, lenNN = len(NN)

    hit_count = miss_count = hit_diff = miss_diff = diff = 0.0
    classtype   = var['classType']
    datatype    = attr[header[ai]][0]
    mmdiff      = attr[header[ai]][3]
    inst_item   = x[inst][ai]

    #--------------------------------------------------------------------------
    if(classtype == 'discrete'):
        for i in range(lenNN):
            NN_item = x[NN[i]][ai]
            if(math.isnan(inst_item) or math.isnan(NN_item)): continue
            if(datatype == 'continuous'):
                calc = abs(inst_item - NN_item)/mmdiff

            if(y[inst] == y[NN[i]]):
                hit_count += 1  # HIT
                if(inst_item != NN_item):
                    if(datatype == 'continuous'):
                        hit_diff -= calc
                    else:
                        hit_diff -= 1
            else:  #MISS
                miss_count += 1
                if(inst_item != NN_item):
                    if(datatype == 'continuous'):
                        miss_diff += calc
                    else:
                        miss_diff += 1

        return hit_diff + miss_diff
    #--------------------------------------------------------------------------
    if(classtype == 'multiclass'):
        class_store = dict()
        missClassPSum = 0

        for each in mcmap:
            if(each != y[inst]):
                class_store[each] = [0,0]
                missClassPSum += mcmap[each]

        for i in range(len(NN)):
            NN_item = x[NN[i]][ai]
            if(math.isnan(inst_item) or math.isnan(NN_item)): continue
            if(datatype == 'continuous'):
                calc = abs(inst_item - NN_item)/mmdiff

            if(y[inst] == y[NN[i]]): #HIT
                hit_count += 1
                if(inst_item != NN_item):
                    if(datatype == 'continuous'):
                        hit_diff -= calc
                    else:
                        hit_diff -= 1
            else:  #MISS
                for missClass in class_store:
                    if(y[NN[i]] == missClass):
                        class_store[missClass][0] += 1
                        if(inst_item != NN_item):
                            if(datatype == 'continuous'):
                                class_store[missClass][1] += calc
                            else:
                                class_store[missClass][1] += 1
        #Corrects for both multiple classes, as well as missing data.
        missSum = 0
        for each in class_store:
            missSum += class_store[each][0]
        missAvg = missSum/float(len(class_store))

        hit_prop = hit_count/float(len(NN))  # correct for missing data
        for each in class_store:
            diff += (mcmap[each]/float(missClassPSum)) * class_store[each][1]

        diff = diff * hit_prop
        miss_prop = missAvg/float(len(NN))
        diff += hit_diff * miss_prop

        return diff
    #--------------------------------------------------------------------------
    if(classtype == 'continuous'):

        for i in range(len(NN)):
            NN_item = x[NN[i]][ai]
            if(math.isnan(inst_item) or math.isnan(NN_item)): continue
            if(datatype == 'continuous'):
                calc = abs(inst_item - NN_item)/mmdiff

            if(abs(y[inst] - y[NN[i]]) < same_class_bound):  #HIT
                hit_count += 1
                if(inst_item != NN_item):
                    if(datatype == 'continuous'):
                        hit_diff -= calc
                    else:
                        hit_diff -= 1
            else:   #MISS
                miss_count += 1
                if(inst_item != NN_item):
                    if(datatype == 'continuous'):
                        miss_diff += calc
                    else:
                        miss_diff += 1

        # Take hit/miss inbalance into account (coming from missing data, or
        # inability to find enough continuous neighbors)

        hit_prop = hit_count/float(len(NN))
        miss_prop = miss_count/float(len(NN))

        return hit_diff * miss_prop + miss_diff * hit_prop
###############################################################################
