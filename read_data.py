# Wed Jul 13 12:41:40 EDT 2016
# NOT PART of PACKAGE (example only)
###############################################################################
# Functions: read_data
###############################################################################
import time as tm
import numpy as np
import sys
###############################################################################
# Read in data file into a numpy array (data) and a header
def np_read_data(fname, options):
    import csv
    
    start = tm.time()
    V = options['verbose']
    md = options['missingdata']
    # determine delimiter ----------------#
    fh = open(fname)
    line = fh.readline().rstrip()
    fh.close()
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(line)
    delim = dialect.delimiter
    #-------------------------------------#
    header = line.split(delim)
    # reading into numpy array
    #dat = pd.read_csv(fname, sep=delim, na_values=md)
    data = np.genfromtxt(fname, missing_values=md, skip_header=1, 
                         dtype=np.double, delimiter=delim)
    if(V):
        print(fname + ": data input elapsed time(sec) = " + str(tm.time() - start))
        sys.stdout.flush()
    return header, data
###############################################################################
