import time
import pandas as pd
import numpy as np
from numpy import random as rnd
import matrixprofile as mpx
import matplotlib.pyplot as plt
import math
from paa import paa
import mass_ts as mts
from plotAMP import plotResult
from scipy.stats import zscore
from decimal import Decimal, getcontext
import pyscamp as scamp
from scipy.signal import find_peaks

%matplotlib widget


def momp(T, m, dd, report = True):

    amp, _ , _ , _ = approxMP(T, m, dsr = dd)
    
    absf, absf_loc, _ = findBSF(amp, m)

    bsf, _ = BSFlocalsearch(T,m, dd, absf_loc)
    
    pruned_ts , final_indices, prunning_perc = prune(T, m, absf, bsf, amp)
    pruned_MP = mpx.compute(pruned_ts,m)['mp']
    _ , min_loc, _ = findBSF(pruned_MP, m)
    momp_out = [final_indices[i] for i in min_loc]


    if report:
        print('Prunning Rate: ', 1- round(prunning_perc,4))

        return momp_out, amp, absf, bsf


def findIndices(amp, m, absf, bsf):
    candidates = []
    initset = np.where((amp >= absf) & (amp <= bsf))[0]
    for ii in initset:
        startidx = max(1, ii - 2)
        endidx = min(len(amp), ii + 2)
        candidates.append(list(range(startidx, endidx)))
    candidates = sum(candidates, [])
    candidates = sorted(list(set(candidates)))
    return candidates


def prune(T, m, absf, bsf, amp):
    idx = findIndices(amp, m, absf, bsf)
    pruned_ts = T[idx]
    prunning_perc = len(pruned_ts) / len(T)
    return pruned_ts, idx, round(prunning_perc,3)


def BSFlocalsearch(ts, m, dsr, absf_loc):
    st = time.time()
    offset = dsr - 1
    width = 2*offset + m
    temp_ts = np.zeros((width *len(absf_loc)))
    for ii, loc in enumerate(absf_loc):
        temp_ts[ii*width : (ii+1)*width] = ts[loc- offset : loc + offset + m]

    mp = mpx.compute(temp_ts,m)['mp']
    bsf = np.min(mp)
    end = time.time()
    bsf_time = end - st

    return bsf, bsf_time
    


def findBSF(mp, m):
    st = time.time()
    mp = np.round(mp, 4)
    min_val = np.min(mp)
    min_indices = np.where(mp == min_val)[0]
    min_locations = [min_indices[0]]
    for ii in range(1, len(min_indices)):
        if min_indices[ii] != min_indices[ii-1] + 1:
            min_locations.append(min_indices[ii])
    end = time.time()
    bsf_time = end - st
    return min_val, min_locations, bsf_time



def approxMP(ts, w, dsr):
    st = time.time()
    if dsr > 1 : 
        ts_dsmp = paa(ts, len(ts)//dsr)
        w_dsmp = np.ceil(w/dsr).astype(int)
        dsmp_mp = mpx.compute(ts_dsmp, w_dsmp)['mp']
        amp = np.sqrt(dsr) * np.repeat(dsmp_mp, dsr)
        end = time.time()
        return amp , dsr, ts_dsmp, end-st
    else:
        end = time.time()
        return mpx.compute(ts,w)['mp'], 1, ts, end-st
    

def genData():
    np.random.seed(0)
    data = rnd.randn(64000)/20
    t = np.arange(0, 10, 0.01)
    motif = np.abs(np.sin(t))
    motif = motif + rnd.randn(len(motif))/20
    print(len(motif))
    data[30000:30000+len(motif)] = motif
    data[42000:42000+len(motif)] = motif
    data = zscore(data + rnd.randn(len(data))/20)
    
    return data


if __name__=="__main__":

    plt.close('all')

    T = genData()
    m = 500
    dd = 8 #Fix down sampling rate for now

    print('T is {} long and the requested motif length is {}'.format(len(T), m))

    st = time.time()
    momp_out, amp, absf, bsf = momp(T, m, dd, report = True)
    end = time.time()
    momp_time = end - st

    st = time.time()
    mp = mpx.compute(T,m)['mp']
    _ , min_loc, _ = findBSF(mp, m)
    mp_time = time.time() - st

    print('Speedup : {}X'.format(round(mp_time / momp_time)))
    print('MOMP: ', momp_out)
    print('MPx: ', min_loc)
    plotResult(T, m, mp, amp, dd, absf, bsf)