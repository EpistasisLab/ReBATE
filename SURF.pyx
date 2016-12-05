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
# this runs both SURF and SURF*
def runSURF(header, double[:,::1] x, double[::1] y, attr, var, distArray, options):
    V = options['verbose']
    algorithm = options['algorithm']
    start = tm.time()
    cdef: 
        int maxInst = var['datalen']
        int numattr = var['NumAttributes']
    Scores = [0] * numattr
    #--------------------------------------------------------------------------
    def find_nearest_neighbor():  # for SURF
        NN = []
        min_indicies = []
        cdef int i

        for i in range(maxInst):
            if(inst != i):
                locator = [inst,i]
                if(i > inst): locator.reverse()
                d = distArray[locator[0]][locator[1]]
                if(d < avgDist):
                    min_indicies.append(i)

        for i in range(len(min_indicies)):
            NN.append(min_indicies[i])

        return NN
    #--------------------------------------------------------------------------
    def find_data_instances():  # for SURFStar
        NN_near=[]; NN_far=[]
        min_indices=[]; max_indices=[]
        cdef int i

        for i in range(maxInst):
            if(inst != i):
                locator = [inst,i]
                if(i > inst): locator.reverse()
                d = distArray[locator[0]][locator[1]]

                if(d < avgDist): min_indices.append(i)
                if(d > avgDist): max_indices.append(i)

        for i in range(len(min_indices)):
            NN_near.append(min_indices[i])

        for i in range(len(max_indices)):
            NN_far.append(max_indices[i])

        return NN_near, NN_far
    #--------------------------------------------------------------------------
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
    #--------------------------------------------------------------------------
    pname = var['phenoTypeName']
    cdef int inst, i

    # calculate avgDist
    alg = options['algorithm']
    sm = cnt = 0
    for i in range(maxInst):
        sm += sum(distArray[i])
        cnt += len(distArray[i])
    avgDist = sm/float(cnt)
    if(V): print('Average Distance = ' + str(avgDist))
    #------------------------------#

    if(var['classType'] == 'multiclass'):
        mcmap = getMultiClassMap()
    else:
        mcmap = 0

    for inst in range(maxInst):
        if(algorithm == 'surf'):
            NN = find_nearest_neighbor()
            NN = np.array(NN, dtype=np.int32)
            if(len(NN) <= 0): continue
            for f in range(var['NumAttributes']):
                Scores[f] += \
                    evaluate_SURF(x,y,header,attr,var,NN,f,inst,mcmap,alg)
        elif(algorithm == 'surfstar'):
            NN_near,NN_far = find_data_instances()
            NN_near = np.array(NN_near, dtype=np.int32)
            NN_far = np.array(NN_far, dtype=np.int32)
            for f in range(var['NumAttributes']):
                if(len(NN_near) > 0):
                    Scores[f] += \
                    evaluate_SURF(x,y,header,attr,var,NN_near,f,inst,mcmap,alg)
                if(len(NN_far) > 0):
                    Scores[f] -= \
                    evaluate_SURF(x,y,header,attr,var,NN_far,f,inst,mcmap,alg)

    if(V): print("surf time = " + str(tm.time() - start))
    return Scores
###############################################################################
# evaluates both SURF and SURF* scores
cdef double evaluate_SURF(double[:,::1] x,double[::1] y, header, attr, var, 
                          int[::1] NN, int feature, int inst, mcmap, algorithm):

    fname = header[feature]
    ftype = attr[fname][0]  # feature type
    ctype = var['classType']  # class type
    cdef:
        double diff_hit = 0
        double diff_miss = 0 
        double count_hit = 0
        double count_miss = 0
        double mmdiff = 1
        double xNNifeature, xinstfeature, diff = 0
        int i

    xinstfeature = x[inst][feature]

    if(ftype == 'continuous'): mmdiff = attr[fname][3]

    if(ctype == 'multiclass'):
        class_Store = dict()
        missClassPSum = 0   # for SURF
        for each in mcmap:
            if(each != y[inst]):
                class_Store[each] = [0,0]
                missClassPSum += mcmap[each]  # for SURF

        for i in range(len(NN)):
            NN[i] = int(NN[i])
            xNNifeature = x[NN[i]][feature]
            absvalue = abs(xinstfeature - xNNifeature)/mmdiff

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

        # corrects for both multiple classes as well as missing data
        missSum = 0
        for each in class_Store: missSum += class_Store[each][0]
        missAverage = missSum/float(len(class_Store))

        hit_proportion = count_hit / float(len(NN)) # Correct for NA
        for each in class_Store:
            if(algorithm == 'surf'):
                diff_miss += (mcmap[each] / float(missClassPSum)) * \
                    class_Store[each][1]
            else:  #surfstar
                diff_miss += (class_Store[each][0]/float(missSum)) * \
                    class_Store[each][1]

        diff = diff_miss * hit_proportion
        miss_proportion = missAverage / float(len(NN))
        diff += diff_hit * miss_proportion

    #--------------------------------------------------------------------------
    elif(ctype == 'discrete'):
        for i in range(len(NN)):
            xNNifeature = x[NN[i]][feature]
            xinstfeature = x[inst][feature]
            absvalue = abs(xinstfeature - xNNifeature)/mmdiff

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

        hit_proportion = count_hit/float(len(NN))
        miss_proportion = count_miss/float(len(NN))
        diff = diff_hit * miss_proportion + diff_miss * hit_proportion
    #--------------------------------------------------------------------------
    else: # CONTINUOUS endpoint
        same_class_bound = var['phenSD']
        for i in range(len(NN)):
            xNNifeature = x[NN[i]][feature]
            xinstfeature = x[inst][feature]
            absvalue = abs(xinstfeature - xNNifeature)/mmdiff

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

        hit_proportion = count_hit/float(len(NN))
        miss_proportion = count_miss/float(len(NN))
        diff = diff_hit * miss_proportion + diff_miss * hit_proportion

    return diff
