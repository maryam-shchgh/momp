
import numpy as np
import matrixprofile as mpx
from paa import paa
import pandas as pd
#For testing purposes
import time 
import matplotlib.pyplot as plt



def approxMP(ts, m, dsr):
    if dsr > 1 : 
        ts_dsmp = paa(ts, len(ts)//dsr)
        w_dsmp = np.ceil(m/dsr).astype(int)
        dsmp_mp = mpx.compute(ts_dsmp, w_dsmp)['mp']
        amp = np.sqrt(dsr) * np.repeat(dsmp_mp, dsr)
        return amp
    else:
        return mpx.compute(ts,m)['mp']

def bsfMotif(mp):

    min_val = np.min(mp)
    min_indices = np.where(mp == min_val)[0]
    min_locations = [min_indices[0]]
    for ii in range(1, len(min_indices)):
        if min_indices[ii] != min_indices[ii-1] + 1:
            min_locations.append(min_indices[ii])

    return min_val, min_locations


def refine(ts, m, dsr, absf_loc, bsf, bsf_loc, idxList):
    offset = dsr - 1
    sub1 = ts[absf_loc[0] : min(absf_loc[0] + m, len(ts))]
    sub2 = ts[absf_loc[1] : min(absf_loc[1] + m, len(ts))]
    sub1 = ts[absf_loc[0] : absf_loc[0] + m]
    sub2 = ts[absf_loc[1] : absf_loc[1] + m]
    bsf_origin = round(np.sqrt(np.sum((sub1-sub2)**2)), 2)

    loc1_st, loc1_end = max(0, absf_loc[0] - offset), min(len(ts), absf_loc[0] + offset)
    loc2_st, loc2_end = max(0, absf_loc[1] - offset), min(len(ts), absf_loc[1] + offset)

    for ii in range(loc1_st , loc1_end):
        sub1 = ts[ii : min(ii + m, len(ts))]
        for jj in range(loc2_st , loc2_end):
            sub2 = ts[jj : min(jj + m, len(ts))]
            temp_bsf = round(np.linalg.norm(sub1-sub2),2)
            if temp_bsf <= bsf:
                bsf = temp_bsf
                bsf_loc = [idxList[ii], idxList[jj]]

    return bsf, bsf_loc, bsf_origin


def prune(T, m, absf, bsf, amp, idxList):
    idx = findIndices(amp, m, absf, bsf, len(T))
    pruned_ts, pruned_idxList = T[idx], idxList[idx]
    return pruned_ts, pruned_idxList


def findIndices(amp, m, absf, bsf, n):

    initset = np.where((amp >= absf) & (amp <= bsf))[0]
    split_indices = np.where(np.diff(initset) != 1)[0] + 1
    subarrays = np.split(initset, split_indices)

    for ii in range(len(subarrays)):
        lastitem = subarrays[ii][-1]
        if ii < len(subarrays) - 1:
            firstitem = subarrays[ii+1][0]
        else:
            firstitem = n
        gap = firstitem  - lastitem
        if gap < m:
            subarrays[ii] = np.concatenate((subarrays[ii], np.arange(lastitem+1, lastitem + gap)))
        else:
            subarrays[ii] = np.concatenate((subarrays[ii], np.arange(lastitem+1, lastitem + m)))

    candidates = np.concatenate(subarrays)

    return candidates

