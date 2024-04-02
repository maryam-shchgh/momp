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

    bsf = np.inf
    momp_out, amp, absf, bsf, _ = process(T, m, dd, np.arange(len(T)), bsf)

    return momp_out, amp, absf, bsf

def process(T, m, dd, idxList, bsf):

    # idxList : original indices for the input data points
    amp, _ , _ , _ = approxMP(T, m, dsr = dd)
    absf, original_bsf_loc= bsfMotif(amp, m)
    bsf,momp_out, original_bsf = exactLocalSearch(T,m, dd, original_bsf_loc, bsf, idxList)
    pruned_ts , idxList, prunning_perc = prune(T, m, absf, bsf, amp, idxList, n_orig)

    print('MOMP: Tpaa1in{} |BSFdist : {} |Loc: {} |'\
          'OriginBSFdist: {} |OriginBSFloc: {} |PrunePerc: {}'\
            .format(dd, round(bsf,2), momp_out,\
                round(original_bsf,2), original_bsf_loc,\
                                        round(prunning_perc,4)))

    if dd//2 >= 4:
       momp_out, amp, absf, bsf, pruned_ts = process(pruned_ts, m, dd//2, idxList, bsf)

    return momp_out, amp, absf, bsf, pruned_ts


def findMaxPow2Sub(N):
 
    # if N is a power of two simply return it
    if (not (N & (N - 1))):
        return N
         
    # else set only the most significant bit
    return 0x8000000000000000 >>  (64 - N.bit_length())

def next_closest_multiple_greater(N, K):
    quotient, remainder = divmod(N, K)
    if remainder == 0:
        return (quotient + 1) * K
    else:
        return (quotient + 1) * K


def main():

    plt.close('all')

    T , m_orig = genData() , 800
    global n_orig
    n_orig = len(T)
    m = findMaxPow2Sub(m_orig)
    dd = m // 4

    print('T length: {} | Requested motif length : {}'.format(n_orig, m_orig))
    print('The largest power-of-two <  {} : {}'.format(m_orig, m))
    print('Starting downsampling rate (1/4 * {}):  {}'.format(m, m//4))

    
    padded_length = next_closest_multiple_greater(n_orig, dd)
    padcount = padded_length - n_orig
    T = np.pad(T, (0, padcount), mode='constant', constant_values=0)

    #Evaluation of the result
    st = time.time()
    mp = mpx.compute(T,m)['mp']
    _, min_loc = bsfMotif(mp, m)
    mp_time = round(time.time() - st,2)
    print('MPx: {} | minval : {}'.format(min_loc, round(np.min(mp),2)))

    st = time.time()
    _, amp, absf, bsf = momp(T, m, dd)
    momp_time = round(time.time() - st,2)

    print('Speedup : {}X | MPx:{}s | MOMP: {}s'.format(round(mp_time / momp_time), mp_time, momp_time))
    plotResult(T, m, mp, amp, absf, bsf, dd)


if __name__=="__main__":

    # Suppress the warning
    warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    main()

    