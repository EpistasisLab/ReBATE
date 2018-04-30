# Fri Jun 24 14:28:10 EDT 2016
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
import sys
import numpy as np
from math import isnan
#from rebateCommon import getMultiClassMap
###############################################################################
# get multiSURFStar scores for multiclass Class
def runMultiSURFStar(header, double[:,::1] x, double[::1] y, attr, var, distArray, 
                 options):
    """ Controls major MultiSURFStar loops. """

    V = options['verbose']
    start = tm.time()
    #--------------------------------------------------------------------------
    def get_individual_distances():
        d=[]
        for j in range(datalen):
            if (i!=j):
                locator = [i,j]
                if(i < j): locator.reverse()
                d.append(distArray[locator[0]][locator[1]])
        return d
    #--------------------------------------------------------------------------
    def makeClassPairMap():
    #finding number of classes in the dataset and storing them into the map
        classPair_map = dict()
        for each in multiclass_map:
            for other in multiclass_map:
                if(each != other):
                    locator = [each,other]
                    if(each < other): locator.reverse()
                    #locator = sorted(locator, reverse = True)
                    tempString = str(locator[0]) + str(locator[1])
                    if (not classPair_map.has_key(tempString)):
                        classPair_map[tempString] = [0,0]
        return classPair_map
    #--------------------------------------------------------------------------
    # Find number of classes in the dataset and store them into the map
    # Called from relieff, surf and surf*
    def getMultiClassMap():
        mcmap = dict()
        maxInst = datalen

        for i in range(maxInst):
            if(y[i] not in mcmap):
                mcmap[y[i]] = 0
            else:
                mcmap[y[i]] += 1

        for each in var['phenoTypeList']:
            mcmap[each] = mcmap[each]/float(maxInst)
            
        return mcmap
    #--------------------------------------------------------------------------
    cdef: 
        int numattr = var['NumAttributes']
        int datalen = var['datalen']
        int i, j, k
        double mmdiff = 0, calc = 0, calc1 = 0, missSum
        double count_hit_near, count_miss_near
        double count_hit_far, count_miss_far
        double diff_hit_near, diff_miss_near
        double diff_hit_far, diff_miss_far
        double diff, hit_proportion, miss_proportion
        double xik, xjk

    ScoreList=[0] * numattr
    D=[]; avg_dist=[]
    
    multiclass_map = getMultiClassMap()

    for i in range(datalen):
        dist_vect = get_individual_distances()
        avg_dist.append(np.average(dist_vect))
        D.append(np.std(dist_vect)/2.0)
        
    if(V):
        print("get individual distances elapsed time = " + str(tm.time() - start))
        sys.stdout.flush()

    start = tm.time()
    for k in range(numattr):    #looping through attributes
        datatype = attr[header[k]][0]
        if(datatype == 'continuous'):
            mmdiff = attr[header[k]][3]
            
        count_hit_near = count_miss_near = 0
        count_hit_far = count_miss_far = 0
        diff_hit_near = diff_miss_near = 0
        diff_hit_far = diff_miss_far = 0
        
        class_Store_near = makeClassPairMap()
        class_Store_far = makeClassPairMap()
        
        for i in range(datalen):                     
            xik = x[i][k]
            if(isnan(xik)): continue
            for j in range(i,datalen):
                xjk = x[j][k]
                if(i == j or isnan(xjk)): continue

                if(datatype == 'continuous'):
                    calc = abs(xik - xjk) / mmdiff
                    calc1 = (1 - abs(xik - xjk)) / mmdiff

                locator = [i,j]
                if(i < j): locator.reverse()
                d = distArray[locator[0]][locator[1]]
                #--------------------------------------------------------------
                if (d < (avg_dist[i] - D[i])): #Near
                    if(y[i] == y[j]):
                        count_hit_near += 1
                        if(xik != xjk):
                            if(datatype == 'continuous'):
                                diff_hit_near -= calc
                            else:
                                diff_hit_near -= 1
                    else:
                        count_miss_near += 1
                        locator = [y[i],y[j]]
                        if(y[i] < y[j]): locator.reverse()
                        tempString = str(locator[0]) + str(locator[1])
                        class_Store_near[tempString][0] += 1
                        if(xik !=xjk):
                            if(datatype == 'continuous'):
                                class_Store_near[tempString][1] += calc
                            else:#Discrete
                                class_Store_near[tempString][1] += 1
                #--------------------------------------------------------------
                if (d > (avg_dist[i] + D[i])): #Far
                        
                    if(y[i] == y[j]):
                        count_hit_far += 1
                        if(datatype == 'continuous'):
                            diff_hit_far -= calc1  # Attribute being similar is
                        else:                      # more important.
                            if(xik == xjk):
                                diff_hit_far -= 1
                    else:
                        count_miss_far += 1
                        locator = [y[i],y[j]]
                        if(y[i] < y[j]): locator.reverse()
                        tempString = str(locator[0]) + str(locator[1])
                        class_Store_far[tempString][0] += 1
                        
                        if(datatype == 'continuous'):
                            class_Store_far[tempString][1] += calc 
                        else:
                            if(xik == xjk):
                                class_Store_far[tempString][1] += 1    
        #Near
        missSum = 0 
        for each in class_Store_near:
            missSum += class_Store_near[each][0]
                         
        hit_proportion = count_hit_near/float(count_hit_near+count_miss_near) 
        miss_proportion = count_miss_near/float(count_hit_near+count_miss_near) 
        
        for each in class_Store_near:
            diff_miss_near += (class_Store_near[each][0]/float(missSum))*class_Store_near[each][1]
        diff_miss_near = diff_miss_near * float(len(class_Store_near))

        diff = diff_miss_near*hit_proportion + diff_hit_near*miss_proportion 
                         
        #Far
        missSum = 0 
        for each in class_Store_far:
            missSum += class_Store_far[each][0]

        hit_proportion = count_hit_far/float(count_hit_far + count_miss_far)
        miss_proportion = count_miss_far/float(count_hit_far + count_miss_far) 
        
        for each in class_Store_far:
            diff_miss_far += (class_Store_far[each][0]/float(missSum))*class_Store_far[each][1]

        diff_miss_far = diff_miss_far * float(len(class_Store_far))
        
        diff += diff_miss_far*hit_proportion + diff_hit_far*miss_proportion   
           
        ScoreList[k] += diff    
                
        
    if(V):
        print("MultiSURFStar scoring elapsed time = " + str(tm.time() - start))
        sys.stdout.flush()
    return ScoreList
