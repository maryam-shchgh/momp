
import time
import numpy as np
import matrixprofile as mpx
from paa import paa


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


def prune(T, absf, bsf, amp, idxList, n_orig):
    idx = np.where((amp >= absf) & (amp <= bsf))[0]
    pruned_ts, idxList = T[idx], idxList[idx]
    remaining_perc = len(pruned_ts) / n_orig
    return pruned_ts, idxList, 1- round(remaining_perc,4)


def exactLocalSearch(ts, m, dsr, absf_loc):
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
    


def bsfMotif(mp, m):
    st = time.time()
    # mp = np.round(mp, 4)
    min_val = np.min(mp)
    min_indices = np.where(mp == min_val)[0]
    min_locations = [min_indices[0]]
    for ii in range(1, len(min_indices)):
        if min_indices[ii] != min_indices[ii-1] + 1:
            min_locations.append(min_indices[ii])
    end = time.time()
    bsf_time = end - st
    return min_val, min_locations, bsf_time



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






