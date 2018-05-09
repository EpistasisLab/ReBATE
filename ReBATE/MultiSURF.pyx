# Sun Jul 10 18:29:28 EDT 2016
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
from math import isnan
###############################################################################
# this runs both MultiSURF and MultiSURF*
def runMultiSURF(header, double[:,::1] x, double[::1] y, attr, var, distArray, options):
    V = options['verbose']
    algorithm = options['algorithm']
    avg_dist = []
    D = []
    start = tm.time()
    cdef: 
        int datalen = var['datalen']
        int numattr = var['NumAttributes']
    Scores = [0] * numattr
    #--------------------------------------------------------------------------
    def get_individual_distances():
        d=[]
        for j in range(datalen):
            if (inst!=j):
                locator = [inst,j]
                if(inst < j): locator.reverse()
                d.append(distArray[locator[0]][locator[1]])
        return d
    #--------------------------------------------------------------------------
    def find_nearest_neighbor():  # for MultiSURF
        NN = []
        cdef int j
        for j in range(datalen):
            if(inst != j):
                locator = [inst,j]
                if(j > inst): locator.reverse()
                d = distArray[locator[0]][locator[1]]
                if (d < (avg_dist[inst] - D[inst])): #Near
                    NN.append(j)
        return NN
    #--------------------------------------------------------------------------
    def find_data_instances():  # for MultiSURFStar
        NN_near=[]; NN_far=[]
        min_indices=[]; max_indices=[]
        cdef int i
        for j in range(datalen):
            if(inst != j):
                locator = [inst,j]
                if(j > inst): locator.reverse()
                d = distArray[locator[0]][locator[1]]
                if (d < (avg_dist[inst] - D[inst])): #Near
                    NN_near.append(j)
                if (d > (avg_dist[inst] + D[inst])): #Far
                    NN_far.append(j)
        return NN_near, NN_far
    #--------------------------------------------------------------------------
    # Find number of classes in the dataset and store them into the map
    def getMultiClassMap():
        mcmap = dict()

        for i in range(datalen):
            if(y[i] not in mcmap):
                mcmap[y[i]] = 0
            else:
                mcmap[y[i]] += 1
            
        return mcmap
    #--------------------------------------------------------------------------
    pname = var['phenoTypeName']
    cdef int inst, i

    alg = options['algorithm']

    if(var['classType'] == 'multiclass'):
        mcmap = getMultiClassMap()
    else:
        mcmap = 0

    for inst in range(datalen):
        #Precalculate average distances and sd's on a 'per target instance' basis (important to MultiSURF functionality)
        dist_vect = get_individual_distances()
        avg_dist.append(np.average(dist_vect))
        D.append(np.std(dist_vect)/2.0)
        
        if(algorithm == 'multisurf'):
            NN = find_nearest_neighbor()
            NN = np.array(NN, dtype=np.int32)
            if(len(NN) <= 0): continue
            for f in range(var['NumAttributes']):
                Scores[f] += evaluate_MultiSURF(x,y,header,attr,var,NN,f,inst,mcmap)
                
        elif(algorithm == 'multisurfstar'):
            NN_near,NN_far = find_data_instances()
            NN_near = np.array(NN_near, dtype=np.int32)
            NN_far = np.array(NN_far, dtype=np.int32)
            for f in range(var['NumAttributes']):
                if(len(NN_near) > 0):
                    Scores[f] += evaluate_MultiSURF(x,y,header,attr,var,NN_near,f,inst,mcmap)
                if(len(NN_far) > 0):
                    Scores[f] += evaluate_MultiSURF_Far(x,y,header,attr,var,NN_far,f,inst,mcmap)

    if(V): print("surf time = " + str(tm.time() - start))
    return Scores
###############################################################################
# evaluates both MultiSURF and MultiSURF* scores for nearest neighbors
cdef double evaluate_MultiSURF(double[:,::1] x,double[::1] y, header, attr, var, 
                          int[::1] NN, int feature, int inst, mcmap):

    fname = header[feature]
    ftype = attr[fname][0]  # feature type
    fstd  = attr[fname][4]
    ctype = var['classType']  # class type
    cdef:
        double diff_hit = 0
        double diff_miss = 0 
        double count_hit = 0
        double count_miss = 0
        double mmdiff = 1
        double xNNifeature, xinstfeature, diff = 0
        int i
        double datalen = var['datalen']

    xinstfeature = x[inst][feature]

    if(ftype == 'continuous'): mmdiff = attr[fname][3]

    if(ctype == 'multiclass'):
        class_Store = dict()
        for each in mcmap:
            if(each != y[inst]):
                class_Store[each] = [0,0]

        for i in range(len(NN)):
            NN[i] = int(NN[i])
            xNNifeature = x[NN[i]][feature]
            absvalue = ramp_function(mmdiff, fstd, xinstfeature, xNNifeature, var) #abs(xinstfeature - xNNifeature)/mmdiff 

            if(isnan(xinstfeature) or isnan(xNNifeature)): continue

            if(y[inst] == y[NN[i]]):  # HIT
                count_hit += 1
                if(xinstfeature != xNNifeature):
                    if(ftype == 'continuous'):
                        diff_hit -= absvalue
                    else:  # discrete
                        diff_hit -= 1

            else:  # MISS
                for missClass in class_Store:
                    if(y[NN[i]] == missClass):
                        class_Store[missClass][0] += 1
                        if(xinstfeature != xNNifeature):
                            if(ftype == 'continuous'):
                                class_Store[missClass][1] += absvalue
                            else:  # discrete
                                class_Store[missClass][1] += 1

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        #Miss component with 'm' normalization
        for each in class_Store:
            count_miss += class_Store[each][0]
            
        if count_hit == 0.0 and count_miss == 0.0:
            return 0.0
        else:
            if count_miss == 0:
                pass
            else: #Normal diff normalization
                for each in class_Store: #multiclass normalization
                    diff += class_Store[each][1] * (class_Store[each][0] / count_miss) * len(class_Store) # Contribution of given miss class weighted by it's observed frequency within NN set.
                diff = diff / count_miss #'m' normalization
            #Hit component: with 'h' normalization
            if count_hit == 0:
                pass
            else:
                diff += (diff_hit / count_hit)
        
        diff = diff / datalen # 'n' normalization
        return diff

    #--------------------------------------------------------------------------
    elif(ctype == 'binary'):
        for i in range(len(NN)):
            xNNifeature = x[NN[i]][feature]
            xinstfeature = x[inst][feature]
            absvalue = ramp_function(mmdiff, fstd, xinstfeature, xNNifeature, var) #abs(xinstfeature - xNNifeature)/mmdiff 

            if(isnan(xinstfeature) or isnan(xNNifeature)): continue

            if(y[inst] == y[NN[i]]):   # HIT
                count_hit += 1
                if(xinstfeature != xNNifeature):
                    if(ftype == 'continuous'):
                        diff_hit -= absvalue
                    else: # discrete
                        diff_hit -= 1
            else: # MISS
                count_miss += 1
                if(xinstfeature != xNNifeature):
                    if(ftype == 'continuous'):
                        diff_miss += absvalue
                    else: # discrete
                        diff_miss += 1

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        if count_hit == 0.0 or count_miss == 0.0: #Special case, avoid division error
            if count_hit == 0.0 and count_miss == 0.0:
                return 0.0
            elif count_hit == 0.0:
                diff =  (diff_miss / count_miss) / datalen
            else: #count_miss == 0.0
                diff =  (diff_hit / count_hit) / datalen
        else: #Normal diff normalization
            diff = ((diff_hit / count_hit) + (diff_miss / count_miss)) / datalen #normalize by k and n - accounts for missing values by dividing by counts rather than k. 
        
    #--------------------------------------------------------------------------
    else: # CONTINUOUS endpoint
        same_class_bound = var['phenSD']
        for i in range(len(NN)):
            xNNifeature = x[NN[i]][feature]
            xinstfeature = x[inst][feature]
            absvalue = ramp_function(mmdiff, fstd, xinstfeature, xNNifeature, var) #abs(xinstfeature - xNNifeature)/mmdiff 

            if(isnan(xinstfeature) or isnan(xNNifeature)): continue

            if(abs(y[inst] - y[NN[i]]) < same_class_bound): # HIT
                count_hit += 1
                if(xinstfeature != xNNifeature):
                    if(ftype == 'continuous'):
                        diff_hit -= absvalue
                    else: # discrete
                        diff_hit -= 1
            else: # MISS
                count_miss += 1
                if(xinstfeature != xNNifeature):
                    if(ftype == 'continuous'):
                        diff_miss += absvalue
                    else: # discrete
                        diff_miss += 1

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        if count_hit == 0.0 or count_miss == 0.0: #Special case, avoid division error
            if count_hit == 0.0 and count_miss == 0.0:
                return 0.0
            elif count_hit == 0.0:
                diff =  (diff_miss / count_miss) / datalen
            else: #count_miss == 0.0
                diff =  (diff_hit / count_hit) / datalen
        else: #Normal diff normalization
            diff = ((diff_hit / count_hit) + (diff_miss / count_miss)) / datalen #normalize by k and n - accounts for missing values by dividing by counts rather than k. 

    return diff

# evaluates MultiSURF* scores for far instances.  A special scoring update variation improved from SURF* for more efficient scoring updates.  Focuses on scoring features with the 'same' value.
cdef double evaluate_MultiSURF_Far(double[:,::1] x,double[::1] y, header, attr, var, 
                          int[::1] NN, int feature, int inst, mcmap):

    fname = header[feature]
    ftype = attr[fname][0]  # feature type
    fstd  = attr[fname][4]
    ctype = var['classType']  # class type
    cdef:
        double diff_hit = 0
        double diff_miss = 0 
        double count_hit = 0
        double count_miss = 0
        double mmdiff = 1
        double xNNifeature, xinstfeature, diff = 0
        int i
        double datalen = var['datalen']

    xinstfeature = x[inst][feature]

    if(ftype == 'continuous'): mmdiff = attr[fname][3]

    if(ctype == 'multiclass'):
        class_Store = dict()
        missClassPSum = 0   # for MultiSURF
        for each in mcmap:
            if(each != y[inst]):
                class_Store[each] = [0,0]
                missClassPSum += mcmap[each]  # for MultiSURF

        for i in range(len(NN)):
            NN[i] = int(NN[i])
            xNNifeature = x[NN[i]][feature]
            absvalue = ramp_function(mmdiff, fstd, xinstfeature, xNNifeature, var) #abs(xinstfeature - xNNifeature)/mmdiff 

            if(isnan(xinstfeature) or isnan(xNNifeature)): continue

            if(y[inst] == y[NN[i]]):  # HIT
                count_hit += 1
                if(ftype == 'continuous'): 
                    diff_hit -= (1-absvalue) #Sameness should yield most negative score
                else:  # discrete
                    if(xinstfeature == xNNifeature):
                        diff_hit -= 1

            else:  # MISS
                for missClass in class_Store:
                    if(y[NN[i]] == missClass):
                        class_Store[missClass][0] += 1
                        if(ftype == 'continuous'):
                            class_Store[missClass][1] += (1-absvalue) #Sameness should yield most negative score
                        else:  # discrete
                            if(xinstfeature == xNNifeature):
                                class_Store[missClass][1] += 1

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        #Miss component with 'm' normalization
        for each in class_Store:
            count_miss += class_Store[each][0]
            
        if count_hit == 0.0 and count_miss == 0.0:
            return 0.0
        else:
            if count_miss == 0:
                pass
            else: #Normal diff normalization
                for each in class_Store: #multiclass normalization
                    diff += class_Store[each][1] * (class_Store[each][0] / count_miss) * len(class_Store) # Contribution of given miss class weighted by it's observed frequency within NN set.
                diff = diff / count_miss #'m' normalization
            #Hit component: with 'h' normalization
            if count_hit == 0:
                pass
            else:
                diff += (diff_hit / count_hit)
        
        diff = diff / datalen # 'n' normalization
        return diff

    #--------------------------------------------------------------------------
    elif(ctype == 'binary'):
        for i in range(len(NN)):
            xNNifeature = x[NN[i]][feature]
            xinstfeature = x[inst][feature]
            absvalue = ramp_function(mmdiff, fstd, xinstfeature, xNNifeature, var) #abs(xinstfeature - xNNifeature)/mmdiff 

            if(isnan(xinstfeature) or isnan(xNNifeature)): continue

            if(y[inst] == y[NN[i]]):   # HIT
                count_hit += 1
                if(ftype == 'continuous'):
                    diff_hit -= (1-absvalue) #Sameness should yield most negative score
                else: # discrete
                    if(xinstfeature == xNNifeature):
                        diff_hit -= 1
            else: # MISS
                count_miss += 1
                if(ftype == 'continuous'):
                    diff_miss += (1-absvalue) #Sameness should yield most negative score
                else: # discrete
                    if(xinstfeature == xNNifeature):
                        diff_miss += 1

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        if count_hit == 0.0 or count_miss == 0.0: #Special case, avoid division error
            if count_hit == 0.0 and count_miss == 0.0:
                return 0.0
            elif count_hit == 0.0:
                diff =  (diff_miss / count_miss) / datalen
            else: #count_miss == 0.0
                diff =  (diff_hit / count_hit) / datalen
        else: #Normal diff normalization
            diff = ((diff_hit / count_hit) + (diff_miss / count_miss)) / datalen #normalize by k and n - accounts for missing values by dividing by counts rather than k. 
        
    #--------------------------------------------------------------------------
    else: # CONTINUOUS endpoint
        same_class_bound = var['phenSD']
        for i in range(len(NN)):
            xNNifeature = x[NN[i]][feature]
            xinstfeature = x[inst][feature]
            absvalue = ramp_function(mmdiff, fstd, xinstfeature, xNNifeature, var) #abs(xinstfeature - xNNifeature)/mmdiff 

            if(isnan(xinstfeature) or isnan(xNNifeature)): continue

            if(abs(y[inst] - y[NN[i]]) < same_class_bound): # HIT
                count_hit += 1
                if(ftype == 'continuous'):
                    diff_hit -= (1-absvalue) #Sameness should yield most negative score
                else: # discrete
                    if(xinstfeature == xNNifeature):
                        diff_hit -= 1
            else: # MISS
                count_miss += 1   
                if(ftype == 'continuous'):
                    diff_miss += (1-absvalue) #Sameness should yield most negative score
                else: # discrete
                    if(xinstfeature == xNNifeature):
                        diff_miss += 1

        #Score Normalizations (in ReliefF there should be no concern with zero neighbors of a certain class) This is a concern for the other methods. 
        if count_hit == 0.0 or count_miss == 0.0: #Special case, avoid division error
            if count_hit == 0.0 and count_miss == 0.0:
                return 0.0
            elif count_hit == 0.0:
                diff =  (diff_miss / count_miss) / datalen
            else: #count_miss == 0.0
                diff =  (diff_hit / count_hit) / datalen
        else: #Normal diff normalization
            diff = ((diff_hit / count_hit) + (diff_miss / count_miss)) / datalen #normalize by k and n - accounts for missing values by dividing by counts rather than k. 

    return diff

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