
from genData import genData
from momp_utils import prune
from momp_utils import approxMP
from momp_utils import bsfMotif
from momp_utils import refine
from momp_utils import findIndices
import matplotlib.pyplot as plt
import numpy as np
from genData import generate_random_walk
from scipy.stats import zscore
from momp import momp
import matrixprofile as mpx

# %matplotlib widget

def test_dataGen():
    T = genData(dtype = 'test')
    plt.figure(num=1 , figsize=(4,2))
    plt.plot(T)
    plt.show()


def test_amp(T, m):
    amp = approxMP(T, m, dsr = 2)
    absf, absf_loc= bsfMotif(amp)
    plt.figure(figsize=(4,2))
    plt.plot(amp)
    plt.axhline(y=absf, color = 'r')
    plt.axvline(x=absf_loc[0], color = 'r')
    plt.axvline(x=absf_loc[1], color = 'r')
    plt.show()


def test_refine(T, m):
    # T , m, dd = genData(dtype = 'test'), 10, 2
    dd = 2
    amp = approxMP(T, m, dsr = dd)
    absf, absf_loc= bsfMotif(amp)
    bsf,bsf_loc, _ = refine(T,m, dd, absf_loc, np.inf, None, np.arange(len(T)))
    fig = plt.figure(figsize=(4,2))
    axs = fig.subplots(2,1, sharex=True)
    axs[0].plot(T)
    axs[0].axvline(x = bsf_loc[0], color = 'g')
    axs[0].axvline(x = bsf_loc[1], color = 'g')
    axs[1].plot(amp)
    axs[1].axhline(y=absf, color = 'r')
    axs[1].axhline(y=bsf, color = 'r')
    plt.show()


def test_prune(T, m):
    dd = 2
    idxList = np.arange(len(T))
    amp = approxMP(T, m, dsr = 2)
    absf, absf_loc= bsfMotif(amp)
    bsf,bsf_loc, _ = refine(T,m, dd, absf_loc, np.inf, None, idxList)
    pruned_ts, pruned_idxList = prune(T, m, absf, bsf, amp, idxList)
    fig = plt.figure(figsize=(4,2))
    axs = fig.subplots(2,1, sharex=True)
    axs[0].plot(T)
    axs[0].axvline(x = bsf_loc[0], color = 'g')
    axs[0].axvline(x = bsf_loc[1], color = 'g')
    axs[1].plot(amp)
    axs[1].axhline(y=absf, color = 'r')
    axs[1].axhline(y=bsf, color = 'r')
    plt.show()
    plt.figure(figsize=(4,3))
    plt.plot(pruned_ts)
    plt.xticks(np.arange(0,len(pruned_ts),20), [str(ii) for ii in pruned_idxList[::20]])
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


def test_absf():
    # amp = approxMP(T, m, dsr = 2)
    np.random.seed(0)
    amp = np.random.rand(30)
    amp[10:12] = 0
    amp[24:26] = 0
    idxList = np.arange(35)
    absf, absf_loc= bsfMotif(amp, idxList)
    plt.figure(num=2 , figsize=(4,2))
    plt.plot(amp)
    for ii in range(len(absf_loc)):
        plt.axvline(x=absf_loc[ii], color = 'r')
    plt.show()
    print(absf)
    print(absf_loc)



def main():
    # units: 1- dataGen, 2- amp, 3- localSearch , 4- pruning , 5- findIndices
    unit = 0
    T = generate_random_walk(2000) /10
    T[300:500]  = np.abs(np.sin(np.arange(0,2,0.01))) /4
    T[1000:1200]  = np.abs(np.sin(np.arange(0,2,0.01)))/ 4
    T = zscore(T + np.random.rand(len(T))/20)
    m = 128
    plt.figure(figsize=(4,3))
    plt.plot(T)

    mp = mpx.compute(T,m)['mp']
    min_val , min_loc = bsfMotif(mp)
    print('MP: {}, {}'.format(round(min_val,2), min_loc))
    plt.figure(figsize=(4,3))
    plt.plot(mp)
    plt.show()

    _ = momp(T, m, verbose = 0)



    print(' ====  Testing units  ====')
    if unit == 1:
        test_dataGen()
    if unit == 2:
        test_amp(T, m)
    if unit == 3:
        test_refine(T, m)
    if unit == 4:
        test_prune(T, m)
    if unit == 5:
        test_findIndices()
    if unit == 6:
        test_absf()
    



if __name__=="__main__":
    main()