import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import zscore
import numpy as np
import os
import matrixprofile as mpx
from momp_utils import bsfMotif
from datetime import datetime


def vis(T_origin, T, m, dd, uamp, absf, bsf, bsf_loc, absf_loc, pruned_T, pruning_perc):
    
    current_date = datetime.now().strftime("%m-%d-%Y")
    figpath = f"./fig/{current_date}/"

    if not(os.path.exists(figpath)):
        os.makedirs(figpath)


    mp = mpx.compute(T,m)['mp']
    # _ , exact_loc = bsfMotif(mp)

    fig1  = plt.figure(figsize=(4,2.5))
    axs = fig1.subplots(2,1,sharex=True)
    axs[0].plot(T, color = 'k', linewidth = 1)
    axs[0].set_title('T (without downsampling)', fontsize = 10)
    axs[1].plot(uamp , color = 'b')
    axs[1].plot(mp, color = 'c')
    axs[1].axhline(y=absf , color = 'r')
    axs[1].set_title(f'UAMP (m = {m} - dsr = {dd})', fontsize = 10)
    for loc in absf_loc : axs[1].axvline(loc, color = 'r', label = loc)
    # for loc in bsf_loc : axs[0].axvline(loc, color = 'g', label=  loc)
    axs[1].axhline(y=bsf , color = 'g')
    # for ax in axs:
    #     ax.legend()

    plt.subplots_adjust(hspace=0.75)
    plt.show()

    fig2  = plt.figure(figsize=(4,1.2))
    axs = fig2.subplots(1,1,sharex=True)
    axs.plot(T_origin, color = 'k', linewidth = 1)
    for loc in bsf_loc : axs.axvline(loc, color = 'g', label=  loc)
    axs.set_title(f'bsfLoc in Original T', fontsize = 10)
    plt.show()


    if dd > 1:
        fig3  = plt.figure(figsize=(4,1.2))
        axs = fig3.subplots(1,1,sharex=True)
        # axs[0].plot(T)
        axs.plot(pruned_T, color = 'r', linewidth = 1)
        axs.set_title(f'Pruned T (pruning {pruning_perc}%)', fontsize = 10)
        plt.show()

    