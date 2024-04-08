
import numpy as np
import matrixprofile as mpx
from paa import paa
import pandas as pd
#For testing purposes
import time 
import matplotlib.pyplot as plt


def approxMP(ts, m, dsr):

    ts_dsmp = paa(ts, len(ts)//dsr)
    w_dsmp = np.ceil(m/dsr).astype(int)
    dsmp_mp = mpx.compute(ts_dsmp, w_dsmp)['mp']
    amp = np.sqrt(dsr) * np.repeat(dsmp_mp, dsr)

    # fig = plt.figure(figsize = (4,3))
    # plt.plot(amp)

    return amp

def bsfMotif(mp):
    min_val = np.min(mp)
    # print('absf mp min val: ', min_val)
    # print('absf mp vals: ', mp)
    indices = np.where(mp == min_val)[0]
    # print('absf equal indices: ', indices)
    start_indices = np.where(np.diff(np.insert(indices, 0, -2)) > 1)[0]
    # motif_loc = idxList[indices[start_indices]]
    motif_loc = indices[start_indices]
    # min_locations = [min_indices[0]]
    # for ii in range(1, len(min_indices)):
    #     if min_indices[ii] != min_indices[ii-1] + 1:
    #         min_locations.append(min_indices[ii])

    return min_val, motif_loc


# def refine(ts, m, dsr, absf_loc, bsf, bsf_loc):
def refine(ts, m, dsr, absf_loc, bsf):
    offset = dsr - 1
    sub1 = ts[absf_loc[0] : absf_loc[0] + m]
    sub2 = ts[absf_loc[1] : absf_loc[1] + m]
    temp_mp = mpx.compute(np.concatenate((sub1, sub2)),m)['mp']
    local_bsf , _ = bsfMotif(temp_mp)

    loc1_st, loc1_end = max(0, absf_loc[0] - offset), min(len(ts) - m, absf_loc[0] + offset)
    loc2_st, loc2_end = max(0, absf_loc[1] - offset), min(len(ts) - m, absf_loc[1] + offset)

    # for ii in range(loc1_st , loc1_end):
    #     sub1 = ts[ii : ii + m]
    #     for jj in range(loc2_st , loc2_end):
    #         sub2 = ts[jj : jj + m]
    #         temp_mp = mpx.compute(np.concatenate((sub1, sub2)),m)['mp']
    #         temp_bsf , _ = bsfMotif(temp_mp)
    #         if temp_bsf <= bsf:
    #             bsf = temp_bsf
    #             bsf_loc = [ii, jj]

    sub1 = ts[loc1_st : loc1_end + m ]
    for jj in range(loc2_st , loc2_end):
        sub2 = ts[jj : jj + m]
        temp_mp = mpx.compute(sub1,m, sub2)['mp']
        temp_bsf , temp_loc = bsfMotif(temp_mp)
        if temp_bsf <= bsf:
            bsf = temp_bsf
            bsf_loc = [temp_loc[0] + loc1_st, jj]

    # temp_mp = mpx.compute(np.concatenate((ts[loc1_st : loc1_end + m ], ts[loc2_st : loc2_end + m])),m)['mp']
    # temp_bsf , temp_bsf_loc = bsfMotif(temp_mp)
    # if temp_bsf <= bsf:
    #             bsf = temp_bsf
    #             bsf_loc = [loc1_st + temp_bsf_loc[0], loc2_st + temp_bsf_loc[1]]

    return bsf, bsf_loc, local_bsf


def prune(T, m, absf, bsf, amp, idxList):
    idx = findIndices(amp, m, absf, bsf, len(T))
    pruned_ts, pruned_idxList = T[idx], idxList[idx]
    return pruned_ts, pruned_idxList


def findIndices(amp, m, absf, bsf, n):
    # print('amp in pruning: ', amp)
    updated_bsf = min(bsf + (0.2 * bsf), max(amp))
    # updated_bsf = bsf
    initset = np.where((amp >= absf) & (amp <= updated_bsf))[0]
    # print(initset)
    split_indices = np.where(np.diff(initset) != 1)[0] + 1
    # print(split_indices)
    subarrays = np.split(initset, split_indices)

    for ii in range(len(subarrays)):
        lastitem = subarrays[ii][-1]
        if ii < len(subarrays) - 1:
            firstitem = subarrays[ii+1][0]
        else:
            firstitem = n
        gap = firstitem  - lastitem
        if gap < m + m //4:
            subarrays[ii] = np.concatenate((subarrays[ii], np.arange(lastitem+1, lastitem + gap)))
        else:
            subarrays[ii] = np.concatenate((subarrays[ii], np.arange(lastitem+1, lastitem + (m+m//4))))


        firstitem = subarrays[ii][0]
        if ii > 0:
            lastitem = subarrays[ii-1][-1]
        else:
            lastitem = 0
        gap = firstitem  - lastitem
        if gap < m//4:
            subarrays[ii] = np.concatenate((np.arange(firstitem - gap, firstitem), subarrays[ii]))
        else:
            subarrays[ii] = np.concatenate((np.arange(firstitem - (m//4), firstitem), subarrays[ii]))

    candidates = np.concatenate(subarrays)
    candidates = np.unique(candidates)

    return candidates

