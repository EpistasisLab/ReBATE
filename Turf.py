# Wed Jul 13 12:08:30 EDT 2016
import sys
import time as tm
###############################################################################
def runTurf(header, x, y, attr, var, distArray, pct, iterations, fun, options):
    from operator import itemgetter
    import rebateCommon as cmn
    import pandas as pd
    import numpy as np
    
    lost = dict()
    save_x = x
    V = options['verbose']
    if(V): print('Under TURF Control...')
        
    #--------------------------------------------------------------------------
    def adjust_variables(var, attr):
        c = d = 0
        for key in attr:
            if attr[key][0] == 'continuous':
                c += 1
            else:
                d += 1
              
        var['dpct'] = (float(d) / (d + c) * 100, d)
        var['cpct'] = (float(c) / (d + c) * 100, c)
    #--------------------------------------------------------------------------
    def create_newdata(x):
        d = []
        cnt = 0
        
        for a in table:
            if(cnt >= keepcnt):
                lost[a[0]] = iteration + 1
                del attr[a[0]]  # added to remove lost attributes
                i = header.index(a[0])
                del header[i]   # added to remove lost attribute name from header
                x = np.delete(x, i, axis=1)  # added to remove attribute from data
            cnt += 1

        d.append(options['phenotypename'])  # don't forget the class

        if(V):
            print('Getting new variables, attributes and distance array')
            sys.stdout.flush()

        var = cmn.getVariables(header, x, y, options)
        adjust_variables(var, attr)
        cmn.overallDataType(attr,var,options)

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
        if(var['mdcnt'] > 0 or var['dataType'] == 'mixed'):
            import mmDistance as md
            dtypes, diffs = cmn.dtypeArray(header,attr,var)
            distArray = md.getDistances(x, var, dtypes, diffs)
        else:
            distArray = cmn.getDistances(x, attr, var)

        if(V):
            print("get distance array elapsed time(sec) = " + str(tm.time()-begin))
            sys.stdout.flush()

        return header, x, attr, var, distArray, lost
    #--------------------------------------------------------------------------
    fullscores = dict()
    print("Total Iterations: " + str(iterations))
    for iteration in range(iterations):
        numattr = var['NumAttributes']
        if(V):
            print ("============================================")
            print ("Iteration:  " + str(iteration+1))
            print ("Attributes: " + str(numattr))
            sys.stdout.flush()

        table = []

        Scores = fun(header,x,y,attr,var,distArray,options)

        for j in range(var['NumAttributes']):
            if(header[j] == var['phenoTypeName']): continue
            table.append((header[j], Scores[j]))
            fullscores[header[j]] = (Scores[j])

        table = sorted(table,key=itemgetter(1), reverse=True)

        if(iteration + 1 < iterations):
            keepcnt = int(numattr - numattr * pct)
            header,x,attr,var,distArray,lost = create_newdata(x)

    return Scores,save_x,var,fullscores,lost,table

###############################################################################
