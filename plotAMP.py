import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import zscore
import numpy as np
import os
import matrixprofile as mpx


def plotResult(ts, w, mp, amp, absf, bsf, bsf_loc, bsf_loc_origin, dd):
    figpath = './fig/'
    if not(os.path.exists(figpath)):
        os.makedirs(figpath)

    fig = plt.figure(figsize=(4, 3))
    axs = fig.subplots(2, 1, sharex=False)
    ts = zscore(ts)
    axs[0].plot(ts, linewidth = 0.5, color = 'b')
    # axs[0].axvline(x = bsf_loc_origin[0], color='r', linewidth = 0.5)
    # axs[0].axvline(x = bsf_loc_origin[1], color='r', linewidth = 0.5)
    # axs[0].axvline(x = bsf_loc[0], color='g', linewidth = 0.5)
    # axs[0].axvline(x = bsf_loc[1], color='g', linewidth = 0.5)
    axs[1].plot(amp, linewidth = 0.5, color ='b')
    axs[1].axhline(y=absf, color='r', linewidth = 0.7, label='absf')
    axs[1].axhline(y=bsf, color='g', linewidth = 0.5, label='bsf')
    axs[1].set_xlabel('n')


    # Set x-axis ticks to maximum and minimum values
    ylabels = ['znorm Ts', 'UAMP']
    yvar = [ts, mp, amp]
    x_min, x_mean, x_max = 0, len(ts)//2, len(ts)
    for ii, ax in enumerate(axs):
        ax.set_xticks([x_min, x_mean, x_max])
        y_min, y_max = np.min(yvar[ii]),np.max(yvar[ii])
        ax.set_yticks([int(y_min), int(y_max)])
        ax.set_ylabel(ylabels[ii])
        ax.legend()

    plt.tight_layout()
    plt.rcParams.update({'font.size': 6})
    fig.savefig(os.path.join(figpath, \
                'aMP_n_{}_w_{}_dsr_{}.svg'.format(len(ts), w, dd)))
    # fig2.savefig(os.path.join(figpath, \
    #             'sanityCheck_aMP_MP_n_{}_w_{}_dsr_{}.svg'.format(len(ts), w, dd)))
    plt.show()
    