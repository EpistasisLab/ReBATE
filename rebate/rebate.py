#!/usr/bin/env python
#  REBATE CLI
# Thu Apr  6 13:15:38 CDT 2017
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
###############################################################################
import time as tm
import sys
import os
import IO as io
import Common as cmn
###############################################################################
prog_start = tm.time()
Scores = fullscores = table = lost = 0
#-----------------------------------------------------------------------------#
#
# get arguments from command line
#
options = io.getArguments()
V = options['verbose']
turfpct = int(options['turfpct'])
algorithm = options['algorithm']
if(algorithm != 'relieff' and algorithm != 'surf' and algorithm != 'surfstar'
                          and algorithm != 'multisurfstar'):
    print("algorithm " + algorithm + " is not available")
    print("Use relieff, surf, surfstar or multisurfstar")
    sys.exit(1)

if(V):
    print("-------------- Python Version --------------")
    print(sys.version)
    print("--------------------------------------------")

#-----------------------------------------------------------------------------#
# read data into header and numpy matrix
#
input_file = options['filename']
if(os.path.exists(input_file)):
    header, data = io.np_read_data(input_file,options)
else:
    print("File " + input_file + " does NOT exist!")
    sys.exit(1)
#-----------------------------------------------------------------------------#
# get x (data) and y (class) into contiguous numpy arrays
# remove phenotype name from header
x, y = io.getxy(header, data, options)
#-----------------------------------------------------------------------------#
# if there is test data, test it for compatibility
if(options['testdata'] != None):
    testdata = options['testdata']
    if(os.path.exists(testdata)):
        theader, tdata = io.test_testdata(header, testdata, options)
    else:
        print("File " + testdata + " does NOT exist!")
        sys.exit(2)
#-----------------------------------------------------------------------------#
# collect variables
#      var is a dictionary of the following variables:
#      NumAttributes, discreteLimit, discretePhenotype, labelMissingData
#      phenSD, phenotypeList, phenoTypeName, phenoTypeLoc
#      dataType  (overall data type of data)
#      classType (datatype of Class or phenotype)
#
var  = cmn.getVariables(header, x, y, options)

#sys.exit(99)
#-----------------------------------------------------------------------------#
# collect attribute information
#    attributes is a dictionary of tuples:
#    attributes['ATTRIBUTE'] = ('continuous/discrete', MAX/None, MIN,None)
#
attr = cmn.getAttributeInfo(header, x, var, options)
# 
# create header list (cheader) with must headers that have continuous data
#
cheader = []
for i in header:
    if attr[i][0] == 'continuous':
        cheader.append(i)  
if(V):
    print("---------------  Parameters  ---------------")
    print("datafile:   " + options['basename'])
    print("datatype:   " + var['dataType'])
    print("attributes: " + str(var['NumAttributes']))

    if(var['dataType'] == 'mixed'):
        print("    continuous: " + str(var['cpct'][1]))
        print("    discrete:   " + str(var['dpct'][1]))
    print("instances:  " + str(var['datalen']))
    print("missing:    " + str(var['mdcnt']))
    print("classtype:  " + var['classType'])
    if(var['classType'] == 'multiclass'):
        yset = var['phenoTypeList']
        print("  classes:  " + str(len(yset)))
    print("classname:  " + var['phenoTypeName'])
    print("algorithm:  " + options['algorithm'])
    print("--------------------------------------------")
    sys.stdout.flush()
#-----------------------------------------------------------------------------#
# create distance array and remove intermediate data
# if missing and/or mixed data use the mixedDistance function
#
begin = tm.time()
diffs, cidx, didx = cmn.dtypeArray(header, attr, var)
if(var['mdcnt'] > 0):
    import mmDistance as md
    distArray = md.getDistances(x[:,cidx], x[:,didx], var, diffs[cidx])
    disttype = "missing"
else:
    distArray = cmn.getDistances(x, attr, var, cidx, didx, cheader)
    disttype = "discrete/continuous/mixed"
if(V):
    ctime = "[" + tm.strftime("%H:%M:%S") + "]"
    print(ctime + " " + disttype + " distance array time(sec) = " 
                + str(tm.time()-begin))
    sys.stdout.flush()

#-----------------------------------------------------------------------------#
# get Scores based on algorithm selected (-a and -t)
#
if(turfpct > 0):  # Use TURF
    import Turf as T
    pct = float(turfpct)/100.0
    iterations = int(1/float(pct))

    if(algorithm == 'relieff'):
        import relieff as R
        fun = R.runReliefF

    if(algorithm == 'multisurfstar'):
        if(var['classType'] == 'multiclass'):
            import mcmss as MS
        else:
            import multisurfstar as MS
        fun = MS.runMultiSURFstar
    if(algorithm == 'surf' or algorithm == 'surfstar'):
        import surf as S
        fun = S.runSURF

    Scores,x,var,fullscores,lost,table = \
       T.runTurf(header,x,y,attr,var,distArray,pct,iterations,fun,options)
    options['algorithm'] = algorithm + "-turf"

elif(algorithm == 'relieff'):
    import relieff as R
    Scores = R.runReliefF(header,x,y,attr,var,distArray,options)

elif(algorithm == 'multisurfstar'):
    if(var['classType'] == 'multiclass'):
        import mcmss as MS
    else:
        import multisurfstar as MS
    Scores = MS.runMultiSURFStar(header, x, y, attr, var, distArray, options)

elif(algorithm == 'surf' or algorithm =='surfstar'):
    import surf as S
    Scores = S.runSURF(header, x, y, attr, var, distArray, options)
#-----------------------------------------------------------------------------#
# create new data files of some number of top scored attributes
ordered_attr = io.createScoresFile(header, var, Scores, options,
                                   prog_start, turfpct, table, lost)

if(options['topattr'] > 0): 
    io.create_subset(header, x, y, options, ordered_attr)
if(options['testdata'] != None): 
    io.create_test_subset(tdata, options, ordered_attr)

if(V): 
    ctime = "[" + tm.strftime("%H:%M:%S") + "]"
    print(ctime + " Overall program time(sec) = " + str(tm.time() - prog_start))
###############################################################################
