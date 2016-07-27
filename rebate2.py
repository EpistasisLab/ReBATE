import numpy as np
class ReBATE():
    def __init__(self, traindata=None, algorithm='relieff', verbose=False, knearest=10, 
                 dlimit=10, missing='NA', outdir=None, turf=0, testdata=None):
        """ Place to go to run reliefF algorithms """
        self.traindata = traindata
        self.algorithm = algorithm
        self.verbose = verbose
        self.knearest = knearest
        self.dlimit = dlimit
        self.missing = missing
        self.outdir = outdir
        self.turf = turf
        self.testdata = testdata
