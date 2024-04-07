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



def momp(T, m , verbose = True):

    global T_origin
    T_origin = T
    dd = m // 16
    bsf, bsf_loc = np.inf, None
    idxList = np.arange(len(T))
    vis = 1
    bfmp = mpx.compute(T,m)['mp']
    _ , exact_loc = bsfMotif(bfmp)


    print('T length: {} | m: {} | starting dsr:{}'.format(len(T), m, dd))
    
    while dd:

        # for ii in exact_loc:
        #     if ii in idxList:
        #         print('Paa1in{} PASSED > {} in idxList.'.format(dd,ii))
        #     else:
        #         print('Paa1in{} FAILED > {} NOT in idxList.'.format(dd,ii))

        # Locating absf val/ loc using uamp
        uamp = approxMP(T, m, dd)
        # print('checking values of exact loc in uamp dd {}'.format(dd))
        # for ii in exact_loc:
        #     jj = np.where(idxList == ii)[0]
        #     print('uamp value for exact loc {} : {}'.format(ii, uamp[jj-5:jj+5]))

        absf, absf_loc = bsfMotif(uamp)
        absf_loc_pruned_version = absf_loc
        absf_loc = idxList[absf_loc]
        # print('absf loc: ', absf_loc

        #plotting (sanity check)
        if vis:
            fig1  = plt.figure(figsize=(4,2))
            axs = fig1.subplots(2,1,sharex=True)
            axs[0].plot(T)
            axs[0].set_title('Raw Input Time Series', fontsize = 12)
            # for loc in absf_loc_pruned_version : axs[0].axvline(loc, color = 'r', label=  loc)
            axs[1].plot(uamp)
            axs[1].axhline(y=absf , color = 'r')
            axs[1].set_title('Upsampled Approximate MP', fontsize = 12)
            for loc in absf_loc_pruned_version : axs[1].axvline(loc, color = 'r', label = loc)

        if dd == 1:
            plt.subplots_adjust(hspace=0.75)
            plt.show()
            print('MOMP : Tpaa1in{} | BSF: {} | localBSF: {} | BSF loc: {} | localBSF loc: {} | Pruning : {}%'.\
              format(dd, round(absf,2), round(absf,2), absf_loc, absf_loc, pruning))
            return absf, absf_loc
        # Refinement
        
        bsf,bsf_loc, local_bsf = refine(T_origin,m, dd, absf_loc, bsf, bsf_loc)

        #plotting (sanity check)
        # print('local bsf: {}, bsf: {}'.format(local_bsf, bsf))

        if vis:
            for loc in bsf_loc : axs[0].axvline(loc, color = 'g', label=  loc)
            axs[1].axhline(y=bsf , color = 'g')
            for ax in axs:
                ax.legend()

            plt.subplots_adjust(hspace=0.75)
            plt.show()

        # print('Pruning boundaries: {} : {}'.format(absf, bsf))

        # Pruning
        pruned_T , pruned_idxList = prune(T, m, absf, bsf, uamp, idxList) 
        # print('remaining indices: ',pruned_idxList)

        if vis:
            fig2  = plt.figure(figsize=(4,2))
            axs = fig2.subplots(1,1,sharex=True)
            # axs[0].plot(T)
            axs.plot(pruned_T, color = 'r')
            axs.set_title('Pruned Time Series', fontsize = 14)
            plt.show()

        pruning = round(1 - len(pruned_T) / len(T_origin), 2)

        print('MOMP : Tpaa1in{} | BSF: {} | localBSF: {} | BSF loc: {} | localBSF loc: {} | Pruning : {}%'.\
              format(dd, round(bsf,2), round(local_bsf,2), bsf_loc, absf_loc, pruning))
        
        T, idxList = pruned_T, pruned_idxList
        dd = dd // 2



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

    