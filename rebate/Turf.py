# Tue Aug 16 13:26:42 EDT 2016
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
import time as tm
import datetime as dt

###############################################################################
def runTurf(header, x, y, attr, var, distArray, pct, iterations, fun, options, cmn):
    from operator import itemgetter
    import numpy as np

    lost = dict() #Dictionary storing the header names and iteration lost of all features filtered out in this round by TuRF
    start = tm.time()
    save_x = x
    V = options['verbose']

    if(V): print('Under TURF Control...')

    #--------------------------------------------------------------------------
    def create_newdata(header, x):
        dlist = []
        cnt = 0

        if(V):
            print('Reducing attributes by ' + str(options['turfpct']) + '%')
            sys.stdout.flush()

        #Go through table with feature sorted by decreasing scores, once we hit keepcnt, we start adding to lost. 
        for a in table:
            if(cnt >= keepcnt):
                lost[a[0]] = iteration + 1
                i = header.index(a[0])
                dlist.append(i) #store position of each feature removed in dlist. 
            cnt += 1
        
        #update header and dataset to reflect removal of lowest scoring features. 
        header = np.delete(header,dlist).tolist() #remove orphans from header
        x = np.delete(x,dlist,axis=1) #remove orphaned attributes from data
        x = np.ascontiguousarray(x, dtype=np.double)

        if(V):
            print('Getting new variables, attributes and distance array')
            sys.stdout.flush()
        
        #Redo data survey (which may save time in downstream distance array calculation (depending on dataset)
        var = cmn.getVariables(header, x, y, options)
        attr = cmn.getAttributeInfo(header, x, var, options)

        cheader = []
        for i in header:
            if attr[i][0] == 'continuous':
                cheader.append(i)  
                
        if(V):
            print("---------------  Parameters  ---------------")
            print("datatype:   " + var['dataType'])
            print("attributes: " + str(var['NumAttributes']))

            if(var['dataType'] == 'mixed'):
                print("    continuous: " + str(var['cpct'][1]))
                print("    discrete:   " + str(var['dpct'][1]))
            if(var['mdcnt'] > 0):
                print("missing:    " + str(var['mdcnt']))
            print("--------------------------------------------")
            sys.stdout.flush()

        begin = tm.time()
        diffs, cidx, didx = cmn.dtypeArray(header, attr, var)
        #Calculate distance array based on present feature types and data missingness.
        if(var['mdcnt'] > 0):
            import mmDistance as md
            distArray = md.getDistances(x[:,cidx], x[:,didx], var, diffs[cidx])
            disttype = "missing"
        else:
            distArray = cmn.getDistances(x, attr, var, cidx, didx, cheader)
            disttype = "discrete/continuous/mixed"

        if(V):
            print(disttype + " distance array elapsed time(sec) = " 
                    + str(tm.time()-begin))
            sys.stdout.flush()

        return header, x, attr, var, distArray, lost
    #--------------------------------------------------------------------------
    
    print("Total Iterations: " + str(iterations))
    #Main TuRF loop--------------------
    for iteration in range(iterations):
        numattr = var['NumAttributes']
        if(V):
            print ("============================================")
            print ("Iteration:  " + str(iteration+1))
            print ("Attributes: " + str(numattr))
            sys.stdout.flush()

        table = []
        
        #Run the selected core Relief-based algorithm scoring method.
        Scores = fun(header,x,y,attr,var,distArray,options)

        if(V):
            print('Building scores table...')
            sys.stdout.flush()

        for j in range(var['NumAttributes']):
            table.append([header[j], Scores[j]])

        table = sorted(table,key=itemgetter(1), reverse=True)

        if(iteration + 1 < iterations):
            #Calculate features to preserve in the next score update.
            keepcnt = int(numattr - numattr * pct)
            if keepcnt == numattr: #Special case (Ensure at least one feature filtered out in an iteration)
                keepcnt -= 1                
                
            #Store data subset.
            header,x,attr,var,distArray,lost = create_newdata(header, x)

    if(V):
        print('Turf finished! Overall time: ' + str(tm.time() - start))
        sys.stdout.flush()
    return Scores,save_x,var,lost,table
###############################################################################
