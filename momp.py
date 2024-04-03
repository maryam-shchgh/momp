import time
import pandas as pd
import numpy as np
from numpy import random as rnd
import matrixprofile as mpx
import matplotlib.pyplot as plt
from paa import paa
import mass_ts as mts
from plotAMP import plotResult
# from scipy.stats import zscore
# from decimal import Decimal, getcontext
# import pyscamp as scamp
# import os
from genData import genData
from momp_utils import approxMP
from momp_utils import bsfMotif
from momp_utils import prune
from momp_utils import exactLocalSearch
import warnings

%matplotlib widget


def momp(T, m, dd):

    bsf, bsf_loc = np.inf,None
    idxList = np.arange(len(T))
    while dd > 0.5:
        T, idxList, bsf, bsf_loc = process(T, m, dd, idxList, bsf, bsf_loc)
        dd = dd // 2
        

def process(T, m, dd, idxList, bsf, bsf_loc):
    # idxList : original indices for the input data points
    amp  = approxMP(T, m, dd)
    absf, absf_loc= bsfMotif(amp)
    if dd > 1:
        bsf,bsf_loc, bsf_origin = exactLocalSearch(T,m, dd, absf_loc, bsf, bsf_loc, idxList)
        bsf_loc_origin = idxList[absf_loc]
        pruned_T , pruned_idxList = prune(T, m, absf, bsf, amp, idxList)
        pruning_perc = round(1- (len(pruned_T) / n_orig), 3)
        print('MOMP: Tpaa1in{} |BSFdist : {} |Loc: {} |'\
            'OriginBSFdist: {} |OriginBSFloc: {} |PrunePerc: {}'\
                .format(dd, bsf, bsf_loc,\
                            bsf_origin, bsf_loc_origin,\
                            pruning_perc))
        return pruned_T, pruned_idxList, bsf, bsf_loc
    else:
        pruning_perc = round(1- (len(T) / n_orig), 3)
        absf, absf_loc = round(absf,2), idxList[absf_loc]
        print('MOMP: Tpaa1in{} |BSFdist : {} |Loc: {} |'\
            'OriginBSFdist: {} |OriginBSFloc: {} |PrunePerc: {}'\
                .format(dd, absf, absf_loc,\
                            absf, absf_loc,\
                            pruning_perc))
        return T, idxList, absf, absf_loc


# def findMaxPow2Sub(N):
 
#     # if N is a power of two simply return it
#     if (not (N & (N - 1))):
#         return N
         
#     # else set only the most significant bit
#     return 0x8000000000000000 >>  (64 - N.bit_length())

def next_closest_multiple_greater(N, K):
    quotient, remainder = divmod(N, K)
    if remainder == 0:
        return (quotient + 1) * K
    else:
        return (quotient + 1) * K


def main():

    plt.close('all')

    T , m = genData(dtype='random') , 256
    global n_orig
    n_orig = len(T)
    #m = findMaxPow2Sub(m_orig)
    dd = m // 4

    print('T length: {} | Requested motif length : {}'.format(n_orig, m))
    # print('The largest power-of-two <  {} : {}'.format(m_orig, m))
    print('Starting downsampling rate (1/4 * {}):  {}'.format(m, m//4))

    
    padded_length = next_closest_multiple_greater(n_orig, dd)
    padcount = padded_length - n_orig
    T = np.pad(T, (0, padcount), mode='constant', constant_values=0)

    #Evaluation of the result
    st = time.time()
    global mp
    mp = mpx.compute(T,m)['mp']
    _, min_loc = bsfMotif(mp)
    mp_time = round(time.time() - st,2)
    print('MPx: {} | minval : {}'.format(min_loc, round(np.min(mp),2)))

    st = time.time()
    momp(T, m, dd)
    momp_time = round(time.time() - st,2)

    print('Speedup : {}X | MPx:{}s | MOMP: {}s'.format(round(mp_time / momp_time), mp_time, momp_time))
    # plotResult(T, m, mp, amp, absf, bsf, dd)


if __name__=="__main__":

    # Suppress the warning
    warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    main()

    