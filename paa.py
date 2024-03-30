import numpy as np
import matplotlib.pyplot as plt

def paa(s, nseg):

    #s : sequence vector (N X 1)
    #nseg: number of PAA segments
    #data : PAA seg (N X 1)

    N = len(s)
    segLen = N // nseg   #10000
    sN = np.reshape(s[:nseg*segLen], (segLen, nseg), order = 'F')
    avg = np.mean(sN, axis = 0)
    data = np.tile(avg, (1, 1))  # expand segments: np.tile(avg, (segLen, 1))
    data = data.flatten(order='F')  # make column

    return data
    