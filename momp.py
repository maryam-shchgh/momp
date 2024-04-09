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
from datetime import datetime
from visualization_utils import vis
# from scipy.stats import zscore
import os

%matplotlib widget



def momp(T, m , verbose = True, plotting = True):

    global T_origin
    T_origin = T
    dd = m // 4
    bsf, bsf_loc = np.inf, None
    idxList = np.arange(len(T))
    st = time.time()

    print('T length: {} | m: {} | starting dsr:{}'.format(len(T), m, dd))
    
    while dd:

        # Locating absf val/ loc using uamp
        uamp = approxMP(T, m, dd)
        absf, absf_loc = bsfMotif(uamp)
        absf_loc_pruned_version = absf_loc
        absf_loc = idxList[absf_loc]

        # fig1  = plt.figure(figsize=(4,2.5))
        # axs = fig1.subplots(2,1,sharex=True)
        # axs[0].plot(T, color = 'k', linewidth = 1)
        # axs[0].set_title('T (without downsampling)', fontsize = 10)
        # axs[1].plot(uamp , color = 'b')
        # axs[1].plot(mpx.compute(T,m)['mp'] , color = 'g')
        # axs[1].axhline(y=absf , color = 'r')
        # axs[1].set_title(f'UAMP (m = {m} - dsr = {dd})', fontsize = 10)
        # for loc in absf_loc_pruned_version : axs[1].axvline(loc, color = 'r', label = loc)
        # print('absf val: ', absf)


        if dd == 1:
            end = time.time()
            momp_time = round(end - st,5)
            if plotting : vis(T_origin, T, m, dd, uamp, absf, absf, absf_loc, absf_loc_pruned_version, None, None)
            print('MOMP : Tpaa1in{} | BSF: {} | localBSF: {} | BSF loc: {} | localBSF loc: {} | Pruning : {}%'.\
              format(dd, round(absf,2), round(absf,2), absf_loc, absf_loc, pruning))
            return absf, absf_loc, momp_time
        
        # Refinement
        # bsf,bsf_loc, local_bsf = refine(T_origin,m, dd, absf_loc, bsf, bsf_loc)
        bsf,bsf_loc, local_bsf = refine(T_origin,m, dd, absf_loc, bsf)
        print('absf:',absf,'| bsf: ', bsf, np.min(uamp), np.max(uamp))
        # Pruning
        pruned_T , pruned_idxList = prune(T_origin, m, absf, bsf, uamp, idxList) 

        pruning = round(1 - (len(pruned_T) / len(T_origin)), 2)

        if plotting : vis(T_origin, T, m, dd, uamp, absf, bsf, bsf_loc, absf_loc_pruned_version, pruned_T, pruning)
        print('MOMP : Tpaa1in{} | BSF: {} | localBSF: {} | BSF loc: {} | localBSF loc: {} | Pruning : {}%'.\
              format(dd, round(bsf,2), round(local_bsf,2), bsf_loc, absf_loc, pruning))
        
        T, idxList = pruned_T, pruned_idxList
        dd = dd // 2



def next_closest_multiple_greater(N, K):
    quotient, remainder = divmod(N, K)
    if remainder == 0:
        return N
    else:
        return (quotient + 1) * K


def runMP(T, m, verbose = 1, plotting = 1):
    st = time.time()
    mp = mpx.compute(T,m)['mp']
    mp_val, mp_loc = bsfMotif(mp)
    mp_time = round(time.time() - st,5)

    if plotting:
        current_date = datetime.now().strftime("%m-%d-%Y")
        figpath = f"./fig/{current_date}/"
        if not(os.path.exists(figpath)):
            os.makedirs(figpath)
        fig  = plt.figure(figsize=(4,3))
        axs = fig.subplots(2,1,sharex=True)
        axs[0].plot(T, color  = 'b', linewidth = 1)
        axs[0].set_title(f'T (n = {len(T_origin)})', fontsize = 10)
        axs[1].plot(mp, color  = 'r', linewidth = 1)
        axs[1].set_title(f'MP (m = {m})', fontsize = 10)
        plt.show()
        fig.savefig(os.path.join(figpath, \
                    'Torigin_{}_Tpaa1in1_step0_MP_m_{}.svg'.format(len(T_origin),m)))

    if verbose:
        print('MPx >> {} | loc : {} | time: {}s'.format(round(mp_val,2), mp_loc, mp_time))

    return mp, mp_time



def main(verbose = True):

    plt.close('all')
    global n_orig
    global mp
    global T_origin

    n = 128
    T , m = genData(n, mlen = 8, dtype='rwalk') , 8
    dd = m // 4

    if verbose:
        print('T length: {} | Requested motif length : {}'.format(n, m))
        print('Starting downsampling rate (1/4 * {}):  {}'.format(m, m//4))

    padded_length = next_closest_multiple_greater(n, dd)
    padcount = padded_length - n
    T = np.pad(T, (0, padcount), mode='constant', constant_values=0)
    T_origin, n_orig = T, len(T)

    mp, mp_time = runMP(T, m, verbose = 1, plotting=0)

    _, _, momp_time = momp(T, m, verbose = 1, plotting = 1)

    if verbose:
        print('Speedup : {}X | MPx:{}s | MOMP: {}s'.format(round(mp_time / momp_time), mp_time, momp_time))


if __name__=="__main__":

    warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
    main()

    