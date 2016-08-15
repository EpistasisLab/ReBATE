# Mon Aug 15 15:01:33 EDT 2016
import sys
import time as tm
import datetime as dt
###############################################################################
def runTurf(header, x, y, attr, var, distArray, pct, iterations, fun, options):
    from operator import itemgetter
    import Common as cmn
    import numpy as np

    lost = dict()
    start = tm.time()
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
    def create_newdata(header, x):
        dlist = []
        cnt = 0
        tmp = 0
        hlist = []

        if(V):
            print('Reducing attributes by ' + str(options['turfpct']) + '%')
            sys.stdout.flush()

        for a in table:
            if(cnt >= keepcnt):
                lost[a[0]] = iteration + 1
                hlist.append(a[0])     # append lost attribe names to hlist
                i = header.index(a[0])
                dlist.append(i)
            cnt += 1

        header = np.delete(header,dlist).tolist() #remove orphans from header
        x = np.delete(x,dlist,axis=1) #remove orphaned attributes from data
        x = np.ascontiguousarray(x, dtype=np.double)

        if(V):
            print('Getting new variables, attributes and distance array')
            sys.stdout.flush()

        var = cmn.getVariables(header, x, y, options)
        attr = cmn.getAttributeInfo(header, x, var, options)
        #adjust_variables(var, attr)
        #cmn.overallDataType(attr,var,options)

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
            diffs, cidx, didx = cmn.dtypeArray(header,attr,var)
            distArray = md.getDistances(x[:,cidx], x[:,didx], var, diffs[cidx])
            disttype = "mixed/missing"
        else:
            distArray = cmn.getDistances(x, attr, var)
            disttype = "discrete/continuous"

        if(V):
            print(disttype + " distance array elapsed time(sec) = " 
                    + str(tm.time()-begin))
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

        if(V):
            print('Building scores table...')
            sys.stdout.flush()

        for j in range(var['NumAttributes']):
            #if(header[j] == var['phenoTypeName']): continue
            table.append([header[j], Scores[j]])
            fullscores[header[j]] = (Scores[j])

        table = sorted(table,key=itemgetter(1), reverse=True)

        if(iteration + 1 < iterations):
            keepcnt = int(numattr - numattr * pct)
            header,x,attr,var,distArray,lost = create_newdata(header, x)

    if(V):
        print('Turf finished! Overall time: ' + str(tm.time() - start))
        sys.stdout.flush()
    return Scores,save_x,var,fullscores,lost,table

###############################################################################
