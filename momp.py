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
import pyscamp as scamp
import os
from genData import genData
from momp_utils import approxMP
from momp_utils import bsfMotif
from momp_utils import prune
from momp_utils import exactLocalSearch

%matplotlib widget


def momp(T, m, dd):

    momp_out, amp, absf, bsf, _, _ = process(T, m, dd, np.arange(len(T)))

    return momp_out, amp, absf, bsf

def process(T, m, dd, idxList):

    # idxList : original indices for the input data points

    amp, _ , _ , _ = approxMP(T, m, dsr = dd)
    
    absf, absf_loc, _ = bsfMotif(amp, m)

    bsf, _ = exactLocalSearch(T,m, dd, absf_loc)
    
    pruned_ts , idxList, prunning_perc = prune(T, absf, bsf, amp, idxList, n_orig)
    pruned_MP = mpx.compute(pruned_ts,m)['mp']
    _ , min_loc, _ = bsfMotif(pruned_MP, m)
    momp_out = [idxList[i] for i in min_loc]

    print('MOMP: Tpaa1in{} | BSF distance : {} | Loc: [{}, {}] | Prunning Perc: {}%'\
          .format(dd, round(bsf,2), momp_out[0], momp_out[1], round(prunning_perc,4)))

    if dd//2 >= 2:
       momp_out, amp, absf, bsf, pruned_ts, prunning_perc = process(pruned_ts, m, dd//2, idxList)

    return momp_out, amp, absf, bsf, pruned_ts, prunning_perc


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
    T , m_orig = genData() , 300
    global n_orig
    n_orig = len(T)
    m = findMaxPow2Sub(m_orig)
    dd = m // 4
    print('T length: {} | Requested motif length : {}'.format(n_orig, m_orig))
    print('The largest power-of-two <  {} : {}.'.format(m_orig, m))
    print('Starting downsampling rate (1/4 * {}):  {}'.format(m, m//4))

    # padd T up to X the next multiple of dd
    padded_length = next_closest_multiple_greater(n_orig, dd)
    padcount = padded_length - n_orig
    T = np.pad(T, (0, padcount), mode='constant', constant_values=0)

    #Evaluation of the result
    st = time.time()
    mp = mpx.compute(T,m)['mp']
    _ , min_loc, _ = bsfMotif(mp, m)
    mp_time = time.time() - st
    print('MPx: {} | minval : {}'.format(min_loc, round(np.min(mp),2)))
    # fig = plt.figure(num = 2, figsize=(4, 3))
    # plt.plot(mp)
    # fig.savefig(os.path.join('./fig/','mp.svg'))

    st = time.time()
    momp_out, amp, absf, bsf = momp(T, m, dd)
    end = time.time()
    momp_time = end - st

    print('Speedup : {}X'.format(round(mp_time / momp_time)))
    plotResult(T, m, mp, amp, absf, bsf, dd)


if __name__=="__main__":

    plt.close('all')
    main()

    