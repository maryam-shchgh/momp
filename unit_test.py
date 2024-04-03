
from genData import genData
from momp_utils import prune
from momp_utils import approxMP
from momp_utils import bsfMotif
from momp_utils import refine
from momp_utils import findIndices
import matplotlib.pyplot as plt
import numpy as np



def test_dataGen():
    T = genData(dtype = 'test')
    plt.figure(num=1 , figsize=(4,2))
    plt.plot(T)
    plt.show()


def test_amp():
    T , m = genData(dtype = 'test'), 6
    amp = approxMP(T, m, dsr = 2)
    absf, absf_loc= bsfMotif(amp)
    plt.figure(num=2 , figsize=(4,2))
    plt.plot(amp)
    plt.axhline(y=absf, color = 'r')
    plt.show()


def test_refine():
    T , m, dd = genData(dtype = 'test'), 10, 2
    amp = approxMP(T, m, dsr = 2)
    absf, absf_loc= bsfMotif(amp)
    bsf,bsf_loc, _ = refine(T,m, dd, absf_loc, np.inf, None, np.arange(len(T)))
    fig = plt.figure(num=3 , figsize=(4,2))
    axs = fig.subplots(2,1, sharex=True)
    axs[0].plot(T)
    axs[0].axvline(x = bsf_loc[0], color = 'g')
    axs[0].axvline(x = bsf_loc[1], color = 'g')
    axs[1].plot(amp)
    axs[1].axhline(y=absf, color = 'r')
    axs[1].axhline(y=bsf, color = 'r')
    plt.show()


def test_prune():
    T , m, dd = genData(dtype = 'test'), 10, 2
    idxList = np.arange(len(T))
    amp = approxMP(T, m, dsr = 2)
    absf, absf_loc= bsfMotif(amp)
    bsf,bsf_loc, _ = refine(T,m, dd, absf_loc, np.inf, None, idxList)
    pruned_ts, pruned_idxList = prune(T, m, absf, bsf, amp, idxList)
    fig = plt.figure(num=4 , figsize=(4,3))
    axs = fig.subplots(3,1, sharex=True)
    axs[0].plot(T)
    axs[0].axvline(x = bsf_loc[0], color = 'g')
    axs[0].axvline(x = bsf_loc[1], color = 'g')
    axs[1].plot(amp)
    axs[1].axhline(y=absf, color = 'r')
    axs[1].axhline(y=bsf, color = 'r')
    axs[2].plot(pruned_ts)
    plt.show()
    print('Pruning >> idxList: ', pruned_idxList)
    

def test_findIndices():
    T , m, dd = genData(dtype = 'test'), 10, 2
    idxList = np.arange(len(T))
    amp = approxMP(T, m, dsr = 2)
    absf, absf_loc= bsfMotif(amp)
    bsf,_, _ = refine(T,m, dd, absf_loc, np.inf, None, idxList)
    idx = findIndices(amp, m, absf, bsf, len(T))
    # print('findIndices >> Indices:', idx)




if __name__=="__main__":
    # units: 1- dataGen, 2- amp, 3- localSearch , 4- pruning , 5- findIndices
    unit = 5

    print(' ====  Testing units  ====')
    if unit == 1:
        test_dataGen()
    if unit == 2:
        test_amp()
    if unit == 3:
        test_refine()
    if unit == 4:
        test_prune()
    if unit == 5:
        test_findIndices()
    