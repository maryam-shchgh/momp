import time
import pandas as pd
import numpy as np
from numpy import random as rnd
import matrixprofile as mpx
import matplotlib.pyplot as plt
from paa import paa
import mass_ts as mts
from plotAMP import plotResult
from genData import genData
from momp_utils import approxMP
from momp_utils import bsfMotif
from momp_utils import prune
from momp_utils import refine
import warnings
from scipy.stats import zscore
import os

%matplotlib widget

def momp(T, m, verbose = True):
    itertime = 0
    bsf, bsf_loc = np.inf, None
    idxList = np.arange(len(T))
    dd = m // 4
    while dd > 0.5:
        st = time.time()
        T_input = T
        T, idxList, amp, absf, bsf, bsf_loc, bsf_origin, bsf_loc_origin = \
            process(T_input, m, dd, idxList, bsf, bsf_loc)
        pruning_perc = round(1- (len(T) / n_orig), 3)
        end = time.time()
        itertime += (end - st)
        if verbose:
            print('MOMP: Tpaa1in{} |BSFdist : {} |Loc: {} |'\
                'OriginBSFdist: {} |OriginBSFloc: {} |PrunePerc: {} |Time: {}'\
                .format(dd, bsf, bsf_loc, bsf_origin, bsf_loc_origin,\
                        pruning_perc, round(end - st,2)))
            # plotResult(T_input, m, mp, amp, absf, bsf, bsf_loc, bsf_loc_origin, dd)
        
        dd = dd // 2

    return round(itertime,2)
        

def process(T, m, dd, idxList, bsf, bsf_loc):
    amp  = approxMP(T, m, dd)
    absf, absf_loc= bsfMotif(amp)
    if dd > 1:
        bsf,bsf_loc, bsf_origin = refine(T,m, dd, absf_loc, bsf, bsf_loc, idxList)
        bsf_loc_origin = idxList[absf_loc]
        pruned_T , pruned_idxList = prune(T, m, absf, bsf, amp, idxList)
        return pruned_T, pruned_idxList, amp, absf, bsf, bsf_loc, bsf_origin, bsf_loc_origin
    else:
        absf, absf_loc = round(absf,2), idxList[absf_loc]
        return T, idxList, amp, absf, absf, absf_loc, absf, absf_loc


def next_closest_multiple_greater(N, K):
    quotient, remainder = divmod(N, K)
    if remainder == 0:
        return (quotient + 1) * K
    else:
        return (quotient + 1) * K





def main(verbose = True):

    plt.close('all')
    global n_orig
    global mp

    T , m = genData(dtype='random') , 256
    n, dd = len(T), m // 4

    if verbose:
        print('T length: {} | Requested motif length : {}'.format(n, m))
        print('Starting downsampling rate (1/4 * {}):  {}'.format(m, m//4))

    padded_length = next_closest_multiple_greater(n, dd)
    padcount = padded_length - n
    T = np.pad(T, (0, padcount), mode='constant', constant_values=0)
    n_orig = len(T)

    #Evaluation of the result
    st = time.time()
    mp = mpx.compute(T,m)['mp']
    _, min_loc = bsfMotif(mp)
    mp_time = round(time.time() - st,2)
    if verbose:
        print('MPx: {} | minval : {}'.format(min_loc, round(np.min(mp),2)))

    # st = time.time()
    momp_time = momp(T, m, verbose)
    # momp_time = round(time.time() - st,2)

    if verbose:
        print('Speedup : {}X | MPx:{}s | MOMP: {}s'.format(round(mp_time / momp_time), mp_time, momp_time))


if __name__=="__main__":

    warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    main()

    