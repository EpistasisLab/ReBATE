import numpy as np
import getrow as gr
###############################################################################
# This creates a distance array for mixed types of data with or without missing
# data
def getDistances(x, var, dtypes, diffs):
    
    distArray = []
    missing = int(var['mdcnt'])
    
    for index in range(var['datalen']):
        row = gr.get_row_dist(x, dtypes, diffs, index, missing)
        row = list(row)
        distArray.append(row)
        
    return distArray
