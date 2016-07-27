# Mon Jul 18 11:44:19 EDT 2016
import time as tm
import sys
import numpy as np
from math import isnan
###############################################################################
# get multiSURF scores
def runMultiSURF(header, double[:,::1] x, double[::1] y, attr, var, distArray, 
                 options):
    
    V = options['verbose']

    start = tm.time()
    #--------------------------------------------------------------------------
    def get_individual_distances():
        d=[]
        for j in range(len(x)):
            if (i!=j):
                locator = [i,j]
                if(i < j): locator.reverse()
                d.append(distArray[locator[0]][locator[1]])
        return d
    #--------------------------------------------------------------------------
    start = tm.time()
    avg_dist = []; D = []

    cdef:
        int i, j=0, k
        int numattr = var['NumAttributes']
        int datalen = var['datalen']
        double same_class_bound = var['phenSD']
        double xik, xjk, mmdiff = 0.0
        double count_hit_near, count_miss_near, diff_hit_near, diff_miss_near
        double count_hit_far, count_miss_far, diff_hit_far, diff_miss_far
        double d, diff

    pname = var['phenoTypeName']
    classtype = var['classType']
    ScoreList = [0] * numattr

    for i in range(datalen):
        dist_vect = get_individual_distances()
        avg_dist.append(np.average(dist_vect))
        D.append(np.std(dist_vect)/2.0)

    if(V):
        print("get individual distances elapsed time = " + str(tm.time() - start))
        sys.stdout.flush()

    start = tm.time()

    for k in range(numattr):
        datatype = attr[header[k]][0]
        if(datatype == 'continuous'):
            mmdiff = attr[header[k]][3]

        count_hit_near = count_miss_near = diff_hit_near = diff_miss_near  = 0.0
        count_hit_far = count_miss_far = diff_hit_far = diff_miss_far   = 0.0

        for i in range(datalen):
            xik = x[i][k]
            if(isnan(xik)): continue
            for j in range(i,datalen):
                if(i == j) : continue 
                xjk = x[j][k]
                if(isnan(xjk)): continue

                if(datatype == 'continuous'):
                    calc = abs(xik - xjk) / mmdiff

                locator = [i,j]
                if(i < j): locator.reverse()
                d = distArray[locator[0]][locator[1]]

                #--------------------------------------------------------------
                if(d < avg_dist[i] - D[i]):  # NEAR
                    
                    if(classtype == 'discrete'):
                        if(y[i] == y[j]): # SAME ENDPOINT
                            count_hit_near += 1
                            if(xik != xjk):
                                if(datatype == 'continuous'):
                                    diff_hit_near -= calc
                                else:
                                    diff_hit_near -= 1
                        else: # DIFFERENT ENDPOINT
                            count_miss_near += 1
                            if(xik != xjk):
                                if(datatype == 'continuous'):
                                    diff_miss_near += calc
                                else:
                                    diff_miss_near += 1

                    else: # CONTINUOUS ENDPOINT
                        if(abs(y[i] - y[j]) < same_class_bound):
                            count_hit_near += 1
                            if(xik != xjk):
                                if(datatype == 'continuous'):
                                    diff_hit_near -= calc
                                else: # DISCRETE
                                    diff_hit_near -= 1
                        else:
                            count_miss_near += 1
                            if(xik != xjk):
                                if(datatype == 'continuous'):
                                    diff_miss_near += calc
                                else: # DISCRETE
                                    diff_miss_near += 1
                #--------------------------------------------------------------
                if(d > avg_dist[i] + D[i]):  # FAR
                    
                    if(classtype == 'discrete'):
                        if(y[i] == y[j]):
                            count_hit_far += 1
                            if(datatype == 'continuous'):
                                diff_hit_far -= calc
                            else: # DISCRETE
                                if(xik == xjk): diff_hit_far -= 1
                        else:
                            count_miss_far += 1
                            if(datatype == 'continuous'):
                                diff_miss_far += calc
                            else: # DISCRETE
                                if(xik == xjk): diff_miss_far += 1

                    else: # CONTINUOUS ENDPOINT
                        if(abs(y[i] - y[j]) < same_class_bound):
                            count_hit_far += 1
                            if(datatype == 'continuous'):
                                diff_hit_far -= calc
                            else: # DISCRETE
                                if(xik == xjk): diff_hit_far -= 1

                        else:
                            count_miss_far += 1
                            if(datatype == 'continuous'):
                                diff_miss_far += calc
                            else: #DISCRETE
                                if(xik == xjk): diff_miss_far += 1
                #--------------------------------------------------------------

        hit_proportion=count_hit_near/(count_hit_near + count_miss_near)
        miss_proportion=count_miss_near/(count_hit_near + count_miss_near)

        #applying weighting scheme to balance the scores
        diff = diff_hit_near * miss_proportion + diff_miss_near * hit_proportion

        hit_proportion = count_hit_far/(count_hit_far + count_miss_far)
        miss_proportion = count_miss_far/(count_hit_far + count_miss_far)

        #applying weighting scheme to balance the scores
        diff += diff_hit_far * miss_proportion + diff_miss_far * hit_proportion

        ScoreList[k]+=diff

    if(V):
        print("MultiSURF scores elapsed time = " + str(tm.time() - start))
        sys.stdout.flush()

    return ScoreList
###############################################################################
