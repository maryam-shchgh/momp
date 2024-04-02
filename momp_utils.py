
import time
import numpy as np
import matrixprofile as mpx
from paa import paa
import matplotlib.pyplot as plt


def approxMP(ts, m, dsr):
    st = time.time()
    if dsr > 1 : 
        ts_dsmp = paa(ts, len(ts)//dsr)
        w_dsmp = np.ceil(m/dsr).astype(int)
        dsmp_mp = mpx.compute(ts_dsmp, w_dsmp)['mp']
        amp = np.sqrt(dsr) * np.repeat(dsmp_mp, dsr)
        end = time.time()
        return amp , dsr, ts_dsmp, end-st
    else:
        end = time.time()
        return mpx.compute(ts,m)['mp'], 1, ts, end-st


def prune(T, m, absf, bsf, amp, idxList, n_orig):
    idx = findIndices(amp, m, absf, bsf)
    pruned_ts, idxList = T[idx], idxList[idx]
    remaining_perc = len(pruned_ts) / n_orig
    return pruned_ts, idxList, 1- round(remaining_perc,4)


def exactLocalSearch(ts, m, dsr, absf_loc, bsf, idxList):

    offset = dsr - 1

    sub1 = ts[absf_loc[0] : min(absf_loc[0] + m, len(ts))]
    sub2 = ts[absf_loc[1] : min(absf_loc[1] + m, len(ts))]
    original_bsf = np.sqrt(np.sum((sub1-sub2)**2))
    bsf_locs = absf_loc

    loc1_st, loc1_end = max(0, absf_loc[0] - offset), min(len(ts), absf_loc[0] + offset)
    loc2_st, loc2_end = max(0, absf_loc[1] - offset), min(len(ts), absf_loc[1] + offset)

    for ii in range(loc1_st , loc1_end):
        sub1 = ts[ii : min(ii + m, len(ts))]
        for jj in range(loc2_st , loc2_end):
            sub2 = ts[jj : min(jj + m, len(ts))]
            temp_bsf = np.linalg.norm(sub1-sub2)
            if temp_bsf <= bsf:
                bsf = temp_bsf
                bsf_locs = [idxList[ii], idxList[jj]]

    return bsf, bsf_locs, original_bsf
    


def bsfMotif(mp, m):

    # mp = np.round(mp, 4)
    min_val = np.min(mp)
    min_indices = np.where(mp == min_val)[0]
    min_locations = [min_indices[0]]
    for ii in range(1, len(min_indices)):
        if min_indices[ii] != min_indices[ii-1] + 1:
            min_locations.append(min_indices[ii])

    return min_val, min_locations



def findIndices(amp, m, absf, bsf):
    initset = np.where((amp >= absf) & (amp <= bsf))[0]
    # candidates = find_ranges(initset, m)
    candidates = np.array([np.arange(initset[ii], initset[ii] + m) \
                            for ii in range(len(initset))])
    candidates = np.concatenate(candidates)
    candidates = np.unique(candidates, axis=0)
    candidates = candidates[candidates < len(amp)]
    return candidates


def find_ranges(arr, m):

    # break_indices = np.where(np.diff(arr) != 1)[0]
    # slices = np.split(arr, break_indices + 1)
    # consecutive_ranges = np.array([np.arange(slice[0], slice[-1] + m) \
    #                                for slice in slices])
    consecutive_ranges = np.array([np.arange(arr[ii], arr[ii] + m) \
                                   for ii in range(len(arr))])

    return consecutive_ranges






