import time
import pandas as pd
import numpy as np
from numpy import random as rnd
import matrixprofile as mpx
import matplotlib.pyplot as plt
from paa import paa
# import mass_ts as mts
from plotAMP import plotResult
from genData import genData
from momp_utils import approxMP
from momp_utils import bsfMotif
from momp_utils import prune
from momp_utils import refine
import warnings
# from scipy.stats import zscore
# import os

# %matplotlib widget

# def momp(T, m, verbose = True):
#     itertime = 0
#     bsf, bsf_loc = np.inf, None
#     idxList = np.arange(len(T))
#     dd = m // 16
#     while True:
#         st = time.time()
#         T_input = T
#         T, idxList, amp, absf, bsf, bsf_loc, bsf_origin, bsf_loc_origin = \
#             process(T_input, m, dd, idxList, bsf, bsf_loc)
#         pruning_perc = round(1- (len(T) / n_orig), 3)
#         end = time.time()
#         itertime += (end - st)
#         if dd == 1:
#             return   round(itertime,2) 
         
#         if verbose:
#             print('MOMP: Tpaa1in{} |BSFdist : {} |Loc: {} |'\
#                 'OriginBSFdist: {} |OriginBSFloc: {} |PrunePerc: {} |Time: {}'\
#                 .format(dd, bsf, bsf_loc, bsf_origin, bsf_loc_origin,\
#                         pruning_perc, round(end - st,2)))
#             # plotResult(T_input, m, mp, amp, absf, bsf, bsf_loc, bsf_loc_origin, dd)
        
#         dd = dd // 2
#         break



def momp(T, m , verbose = True):

    global T_origin
    T_origin = T
    dd = m // 16
    bsf, bsf_loc = np.inf, None
    idxList = np.arange(len(T))


    print('T length: {} | m: {} | starting dsr:{}'.format(len(T), m, dd))
    
    while dd:

        # Locating absf val/ loc using uamp
        uamp = approxMP(T, m, dd)
        absf, absf_loc = bsfMotif(uamp)
        absf_loc_pruned_version = absf_loc
        absf_loc = idxList[absf_loc]
        # print('absf loc: ', absf_loc)

        print('current idxList start - end :', idxList[0],'-',idxList[-1], '| len : ', len(idxList))
        toprint = [{x: uamp[np.where(idxList == x)[0]]} for x in idxList if ((x > 290 and x < 310) or (x > 990 and x < 1010))] 
        print('Boundary loc: ',toprint)
        if 300 in idxList and 1001 in idxList : print('both motifs are in the boundaries')
        else:
            print('boundary error')

        #plotting (sanity check)
        fig1  = plt.figure(figsize=(4,3))
        axs = fig1.subplots(2,1,sharex=True)
        axs[0].plot(T)
        for loc in absf_loc_pruned_version : axs[0].axvline(loc, color = 'r', label=  loc)
        axs[1].plot(uamp)
        axs[1].axhline(y=absf , color = 'r')
        for loc in absf_loc_pruned_version : axs[1].axvline(loc, color = 'r', label = loc)
        if dd == 1:
            print('MOMP : Tpaa1in{} | BSF: {} | localBSF: {} | BSF loc: {} | localBSF loc: {} | Pruning : {}%'.\
              format(dd, round(absf,2), round(absf,2), absf_loc, absf_loc, pruning))
            return absf, absf_loc
        # Refinement
        
        bsf,bsf_loc, local_bsf = refine(T_origin,m, dd, absf_loc, bsf, bsf_loc)

        #plotting (sanity check)
        # print('local bsf: {}, bsf: {}'.format(local_bsf, bsf))
        for loc in bsf_loc : axs[0].axvline(loc, color = 'g', label=  loc)
        axs[1].axhline(y=bsf , color = 'g')
        for ax in axs:
            ax.legend()

        # Pruning
        pruned_T , pruned_idxList = prune(T, m, absf, bsf, uamp, idxList) 
        # print('remaining indices: ',pruned_idxList)
        fig2  = plt.figure(figsize=(4,2))
        axs = fig2.subplots(1,1,sharex=True)
        # axs[0].plot(T)
        axs.plot(pruned_T, color = 'r')
        plt.show()

        pruning = round(1 - len(pruned_T) / len(T_origin), 2)

        print('MOMP : Tpaa1in{} | BSF: {} | localBSF: {} | BSF loc: {} | localBSF loc: {} | Pruning : {}%'.\
              format(dd, round(bsf,2), round(local_bsf,2), bsf_loc, absf_loc, pruning))
        
        T, idxList = pruned_T, pruned_idxList
        dd = dd // 2

        


# def process(T, m, dd, idxList, bsf, bsf_loc):
#     amp  = approxMP(T, m, dd)
#     absf, absf_loc= bsfMotif(amp)
#     absf_loc = idxList[absf_loc]
#     if dd > 1:
#         # bsf,bsf_loc, bsf_origin = refine(T,m, dd, absf_loc, bsf, bsf_loc, idxList)
#         bsf,bsf_loc, bsf_origin = refine(T_origin,m, dd, absf_loc, bsf, bsf_loc, idxList)
#         if absf != bsf:
#             pruned_T , pruned_idxList = prune(T, m, absf, bsf, amp, idxList)
#         else:
#             pruned_T, pruned_idxList = T, idxList
#         return pruned_T, pruned_idxList, amp, absf, bsf, bsf_loc, bsf_origin, absf_loc
#     else:
#         # absf, absf_loc = round(absf,2), idxList[absf_loc]
#         absf = round(absf,2)
#         return T, idxList, amp, absf, absf, absf_loc, absf, absf_loc


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
    global T_origin

    n = 100000
    T , m = genData(n, dtype='rwalk') , 256
    dd = m // 4

    if verbose:
        print('T length: {} | Requested motif length : {}'.format(n, m))
        print('Starting downsampling rate (1/4 * {}):  {}'.format(m, m//4))

    padded_length = next_closest_multiple_greater(n, dd)
    padcount = padded_length - n
    T = np.pad(T, (0, padcount), mode='constant', constant_values=0)
    T_origin, n_orig = T, len(T)

    #Evaluation of the result
    # st = time.time()
    # mp = mpx.compute(T,m)['mp']
    # _, min_loc = bsfMotif(mp)
    # mp_time = round(time.time() - st,2)

    # fig = plt.figure(figsize=(4,2))
    # axs = fig.subplots(2,1,sharex=True)
    # axs[0].plot(T)
    # axs[1].plot(mp)
    if verbose:
        print('MPx: {} | minval : {}'.format(min_loc, round(np.min(mp),2)))

    momp_time = momp(T, m, verbose)

    if verbose:
        print('Speedup : {}X | MPx:{}s | MOMP: {}s'.format(round(mp_time / momp_time), mp_time, momp_time))


if __name__=="__main__":

    warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    main()

    