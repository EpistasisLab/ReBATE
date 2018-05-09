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
    #--------------------------------------------------------------------------
    def find_nearest_neighbor(indicies):  # for ReliefF
        """Finds the nearest neighbors, but does so differently based on whether 
        endpoint is binary, multiclass, or continuous.  This is only the case for
        ReliefF.  The other implemented ReBATE methods more universally determine
        nearest neighbors. """
        cdef int i, index = 0
        cdef double same_class_bound = var['phenSD']
        NN = []
        if var['classType'] == 'binary' or var['classType'] == 'multiclass':
            for sample in indicies.values(): #sample is all indices with a given class value. 
                n = sorted(sample,key=getsortTuple) #returns sorted list (increasing order)
                NN.extend(n[:neighbors]) #puts in the neighbors up only, not entire sorted list
        else: #continuous endpoint
            #Have to identify hits or misses with respect to differences from target instance's endpoint value.
            #Due to relative nature of 'class' determination the 'class' dictionary must be rebuilt for each target instance. 
            indicies = dict()
            for i in range(2):
                indicies[i] = []
            for i in y: #for each target instance
                if abs(inst-i) < same_class_bound: # considered a hit
                    indicies[0].append(index)
                    index += 1
                else: #considered a miss
                    indicies[1].append(index)
                    index += 1
                    
            for sample in indicies.values(): #sample is all indices with a given class value. 
                n = sorted(sample,key=getsortTuple) #returns sorted list (increasing order)
                NN.extend(n[:neighbors]) #puts in the neighbors up only, not entire sorted list
    
        return NN
    
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
            
        return mcmap
    #-------------------------------------------------------------------------

    cdef:
        int neighbors = var['numNeighbors']
        int numattr = var['NumAttributes']
        int datalen = var['datalen']
        int maxInst = datalen
        
    Scores = [0] * numattr
    
    #Prepare for nearest neighbor identification
    if(var['classType'] == 'multiclass'):
        mcmap = getMultiClassMap()
        indicies = buildIndexLists()
    elif(var['classType'] == 'binary'):
        indicies = buildIndexLists()
        mcmap = 0
    else:
        mcmap = 0
        indicies = None

    #Main Scoring Loop
    for inst in range(datalen): #For each instance as the 'target'
        #Find neighbors
        NN = find_nearest_neighbor(indicies)
        #Update scores for each feature
        for ai in range(numattr):
            #nn = np.array(NN)
            nn = np.ascontiguousarray(NN, dtype=np.int32)
            Scores[ai] += getReliefFScores(header,attr,var,x,y,nn,inst,ai,mcmap)

    #averaging the scores
#     divisor = maxInst * neighbors
#     for ai in range(numattr):
#         Scores[ai] = Scores[ai]/float(divisor)

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
        double datalen = var['datalen']

    hit_count = miss_count = hit_diff = miss_diff = diff = 0.0
    classtype   = var['classType']
    datatype    = attr[header[ai]][0]
    mmdiff      = attr[header[ai]][3]
    fstd        = attr[header[ai]][4]
    inst_item   = x[inst][ai]

    #--------------------------------------------------------------------------
    if(classtype == 'binary'):
        for i in range(lenNN):
            NN_item = x[NN[i]][ai]
            if(math.isnan(inst_item) or math.isnan(NN_item)): continue
            if(datatype == 'continuous'):
                calc = ramp_function(mmdiff, fstd, inst_item, NN_item, var) 

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
                        
        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        if hit_count == 0.0 or miss_count == 0.0: #Special case, avoid division error
            if hit_count == 0.0 and miss_count == 0.0:
                return 0.0
            elif hit_count == 0.0:
                diff =  (miss_diff / miss_count) / datalen
            else: #count_miss == 0.0
                diff =  (hit_diff / hit_count) / datalen
        else: #Normal diff normalization
            diff = ((hit_diff / hit_count) + (miss_diff / miss_count)) / datalen #normalize by k and n - accounts for missing values by dividing by counts rather than k. 
        
        return diff
    #--------------------------------------------------------------------------
    if(classtype == 'multiclass'):
        class_store = dict()

        for each in mcmap:
            if(each != y[inst]):
                class_store[each] = [0,0]

        for i in range(len(NN)):
            NN_item = x[NN[i]][ai]
            if(math.isnan(inst_item) or math.isnan(NN_item)): continue
            if(datatype == 'continuous'):
                calc = ramp_function(mmdiff, fstd, inst_item, NN_item, var) 

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
                                
        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        #Miss component with 'm' normalization
        for each in class_store:
            miss_count += class_store[each][0]
            
        if hit_count == 0.0 and miss_count == 0.0:
            return 0.0
        else:
            if miss_count == 0:
                pass
            else: #Normal diff normalization
                for each in class_store: #multiclass normalization
                    diff += class_store[each][1] * (class_store[each][0] / miss_count) * len(class_store) # Contribution of given miss class weighted by it's observed frequency within NN set.
                diff = diff / miss_count #'m' normalization
            #Hit component: with 'h' normalization
            if hit_count == 0:
                pass
            else:
                diff += (hit_diff / hit_count)
        
        diff = diff / datalen # 'n' normalization
        return diff
    
    #--------------------------------------------------------------------------
    if(classtype == 'continuous'):

        for i in range(len(NN)):
            NN_item = x[NN[i]][ai]
            if(math.isnan(inst_item) or math.isnan(NN_item)): continue
            if(datatype == 'continuous'):
                calc = ramp_function(mmdiff, fstd, inst_item, NN_item, var) 

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

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        if hit_count == 0.0 or miss_count == 0.0: #Special case, avoid division error
            if hit_count == 0.0 and miss_count == 0.0:
                return 0.0
            elif hit_count == 0.0:
                diff =  (miss_diff / miss_count) / datalen
            else: #count_miss == 0.0
                diff =  (hit_diff / hit_count) / datalen
        else: #Normal diff normalization
            diff = ((hit_diff / hit_count) + (miss_diff / miss_count)) / datalen #normalize by k and n - accounts for missing values by dividing by counts rather than k. 
        
        return diff
###############################################################################

def ramp_function(mmdiff, fstd, inst_item, NN_item, var):
    """ Our own user simplified variation of the ramp function suggested by Hong 1994, 1997. Hong's method requires the user to specifiy two thresholds
    that indicate the max difference before a score of 1 is given, as well a min difference before a score of 0 is given, and any in the middle get a 
    score that is the normalized difference between the two continuous feature values. This was done because when discrete and continuous features were mixed,
    continuous feature scores were underestimated.  Towards simplicity, automation, and a dataset adaptable approach, 
    here we simply check whether the difference is greater than the standard deviation for the given feature; if so we assign a score of 1, otherwise we 
    assign the normalized feature score difference.  This should help compensate for the underestimation. """
    #abs(inst_item - NN_item)/mmdiff
    diff = 0
    rawfd = abs(inst_item - NN_item) #prenormalized feature value difference
    
    if var['dataType'] == 'mixed': #Ramp function utilized
        #Check whether feature value difference is greater than the standard deviation
        if rawfd > fstd: #feature value difference is is wider than a standard deviation
            diff = 1
        else:
            diff = rawfd / mmdiff 
        
    else: #Normal continuous feature scoring        
        diff = rawfd / mmdiff 
        
    return diff
#-------------------------------------------------------------------------