import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import zscore
import numpy as np
import os


def plotResult(ts, w, mp, amp, absf, bsf, dd):
    figpath = './fig/'
    if not(os.path.exists(figpath)):
        os.makedirs(figpath)

    fig = plt.figure(num = 1, figsize=(4, 3))
    axs = fig.subplots(2, 1, sharex=True)
    ts = zscore(ts)
    axs[0].plot(ts, linewidth = 0.5)
    axs[1].plot(mp, label='MP', linewidth = 0.5, color = 'b')
    axs[1].plot(amp, label='AMP', linewidth = 0.5, color ='g')
    axs[1].axhline(y=absf, color='r', linewidth = 0.5)
    axs[1].axhline(y=bsf, color='g', linewidth = 0.5)
    axs[1].set_xlabel('n')


    # Set x-axis ticks to maximum and minimum values
    ylabels = ['znorm Ts', 'MP', 'UAMP']
    yvar = [ts, mp, amp]
    x_min, x_mean, x_max = 0, len(ts)//2, len(ts)
    for ii, ax in enumerate(axs):
        ax.set_xticks([x_min, x_mean, x_max])
        y_min, y_max = np.min(yvar[ii]),np.max(yvar[ii])
        ax.set_yticks([int(y_min), int(y_max)])
        ax.set_ylabel(ylabels[ii])
    # for ax in axs.flat:
    #   ax.legend(loc = 'lower left')


    plt.tight_layout()
    plt.rcParams.update({'font.size': 4})
    fig.savefig(os.path.join(figpath, \
                'aMP_n_{}_w_{}_startdsr_{}.svg'.format(len(ts), w, dd)))
    plt.show()
    