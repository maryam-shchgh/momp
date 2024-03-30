import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import zscore
import numpy as np


def plotResult(ts, w, mp, amp, dsr, absf, bsf):
    fig = plt.figure(num = 1, figsize=(4, 3))
    axs = fig.subplots(2, 1, sharex=True)
    fig.suptitle('DSR: {}'.format(dsr))
    ts = zscore(ts)
    axs[0].plot(ts)
    axs[1].plot(mp, label='MP', linewidth = 1, color = 'b')
    axs[1].plot(amp, label='AMP', linewidth = 1, color ='g')
    axs[1].axhline(y=absf, color='r', linewidth = 1)
    axs[1].axhline(y=bsf, color='g', linewidth = 1)
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
    fig.savefig('aMP_n_{}_w_{}_dsr_{}.svg'.format(len(ts), w, dsr))
    plt.show()
    