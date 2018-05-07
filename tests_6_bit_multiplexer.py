
"""
ReBATE was primarily developed at the University of Pennsylvania by:
    - Pete Schmitt (pschmitt@upenn.edu)
    - Ryan J. Urbanowicz (ryanurb@upenn.edu)
    - Weixuan Fu (weixuanf@upenn.edu)
    - and many more generous open source contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

#Initialize hardcoded argument version of rebate.py
import IO as io
import Common as cmn
import relieff as R
import surf as S
import multisurf as MS
import time as tm
import sys
import os

#import warnings

#warnings.filterwarnings('ignore')
###############################################################################

#Setup Options ################################################
options = dict()
options['filename'] = 'data/6Multiplexer_Data_500_0.txt'
options['basename'] = '6Multiplexer_Data_500_0.txt'
options['dir_path'] = 'data'
options['testdata'] = None
options['phenotypename'] = "Class"
options['discretelimit'] = 10
options['neighbors'] = 10
options['missingdata'] = 'NA'
options['algorithm'] = 'relieff'
options['turfpct'] = '0'
options['verbose'] = False
options['debug'] = True
options['topattr'] = 0
options['outputdir'] = '.'
#########################################

#Below is just a copy of the required code from rebate.py#########################
V = options['verbose']
turfpct = int(options['turfpct'])
algorithm = options['algorithm']
if(algorithm != 'relieff' and algorithm != 'surf' and algorithm != 'surfstar' and algorithm != 'multisurfstar' and algorithm != 'multisurf'):
    print("algorithm " + algorithm + " is not available")
    print("Use relieff, surf, surfstar, multisurfstar, or multisurf")
    sys.exit(1)

if(V):
    print("-------------- Python Version --------------")
    print(sys.version)
    print("--------------------------------------------")

#-----------------------------------------------------------------------------#
input_file = options['filename']
if(os.path.exists(input_file)):
    header, data = io.np_read_data(input_file,options)
else:
    print("File " + input_file + " does NOT exist!")
    sys.exit(1)
#-----------------------------------------------------------------------------#
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
var  = cmn.getVariables(header, x, y, options)
attr = cmn.getAttributeInfo(header, x, var, options)
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
#############################################################################
    

###################################################################################################################################################
def test_relieff_Multiplexer():
    """ Test ReliefF on 6-bit Multiplexer"""
    Scores = R.runReliefF(header,x,y,attr,var,distArray,options)
    print("ReliefF + 6-bit Multiplexer ")
    print(str(Scores))
    #Check that score list is not empty
    assert Scores != None
    #Check that a score for all features is output
    assert len(Scores) == 6 #6-bit Multiplexer problem
    #Check that all scores fall between -1 and 1
    assert max(Scores) <= 1 and min(Scores) >= -1    
    #Check that the address bits (indexed as features 0 and 1) have the top scores as expected. 
    indexTopScore = Scores.index(max(Scores))
    assert  indexTopScore == 0 or indexTopScore == 1
    Scores.pop(indexTopScore)
    indexTopScore = Scores.index(max(Scores))
    assert indexTopScore == 0

def test_surf_Multiplexer():
    """ Test SURF on 6-bit Multiplexer"""
    #New parameters
    options['algorithm'] = 'surf'
    Scores = S.runSURF(header, x, y, attr, var, distArray, options)
    print("SURF + 6-bit Multiplexer ")
    print(str(Scores))
    #Check that score list is not empty
    assert Scores != None
    #Check that a score for all features is output
    assert len(Scores) == 6 #6-bit Multiplexer problem
    #Check that all scores fall between -1 and 1
    assert max(Scores) <= 1 and min(Scores) >= -1 
    #Check that the address bits (indexed as features 0 and 1) have the top scores as expected. 
    indexTopScore = Scores.index(max(Scores))
    assert  indexTopScore == 0 or indexTopScore == 1
    Scores.pop(indexTopScore)
    indexTopScore = Scores.index(max(Scores))
    assert indexTopScore == 0
     
def test_surfstar_Multiplexer():
    """ Test SURF* on 6-bit Multiplexer """
    #New parameters
    options['algorithm'] = 'surfstar'
    Scores = S.runSURF(header, x, y, attr, var, distArray, options)
    print("SURF* + 6-bit Multiplexer ")
    print(str(Scores))
    #Check that score list is not empty
    assert Scores != None
    #Check that a score for all features is output
    assert len(Scores) == 6 #6-bit Multiplexer problem
    #Check that all scores fall between -1 and 1
    assert max(Scores) <= 1 and min(Scores) >= -1 
    #Check that the address bits (indexed as features 0 and 1) have the top scores as expected. 
    indexTopScore = Scores.index(max(Scores))
    assert  indexTopScore == 0 or indexTopScore == 1
    Scores.pop(indexTopScore)
    indexTopScore = Scores.index(max(Scores))
    assert indexTopScore == 0
     
def test_multisurfstar_Multiplexer():
    """ Test MultiSURF* on 6-bit Multiplexer """
    #New parameters
    options['algorithm'] = 'multisurfstar'
    Scores = MS.runMultiSURF(header, x, y, attr, var, distArray, options)
    print("MultiSURF* + 6-bit Multiplexer ")
    print(str(Scores))
    #Check that score list is not empty
    assert Scores != None
    #Check that a score for all features is output
    assert len(Scores) == 6 #6-bit Multiplexer problem
    #Check that all scores fall between -1 and 1
    assert max(Scores) <= 1 and min(Scores) >= -1 
    #Check that the address bits (indexed as features 0 and 1) have the top scores as expected. 
    indexTopScore = Scores.index(max(Scores))
    assert  indexTopScore == 0 or indexTopScore == 1
    Scores.pop(indexTopScore)
    indexTopScore = Scores.index(max(Scores))
    assert indexTopScore == 0
    
def test_multisurf_Multiplexer():
    """ Test MultiSURF on 6-bit Multiplexer """
    #New parameters
    options['algorithm'] = 'multisurf'
    Scores = MS.runMultiSURF(header, x, y, attr, var, distArray, options)
    print("MultiSURF + 6-bit Multiplexer ")
    print(str(Scores))
    #Check that score list is not empty
    assert Scores != None
    #Check that a score for all features is output
    assert len(Scores) == 6 #6-bit Multiplexer problem
    #Check that all scores fall between -1 and 1
    assert max(Scores) <= 1 and min(Scores) >= -1 
    #Check that the address bits (indexed as features 0 and 1) have the top scores as expected. 
    indexTopScore = Scores.index(max(Scores))
    assert  indexTopScore == 0 or indexTopScore == 1
    Scores.pop(indexTopScore)
    indexTopScore = Scores.index(max(Scores))
    assert indexTopScore == 0

##################################################################################################################
