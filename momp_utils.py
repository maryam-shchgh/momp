
import numpy as np
import matrixprofile as mpx
from paa import paa
import pandas as pd
#For testing purposes
import time 
import matplotlib.pyplot as plt
from genData import genData


def approxMP(T, m, dd):

    T_dsmp = paa(T, len(T)//dd)
    # m_dsmp = np.ceil(m/dsr).astype(int)
    m_dsmp = m//dd
    # amp = mpx.compute(T_dsmp, m_dsmp)['mp']
    amp = mpx.compute(T_dsmp, m_dsmp)['mp']
    # if dd > 8:
    uamp = np.sqrt(dd) * np.repeat(amp, dd)
    # else:
    #     uamp = np.repeat(amp, dd)
    

    return uamp

def bsfMotif(mp):
    min_val = np.min(mp)
    indices = np.where(mp == min_val)[0]
    start_indices = np.where(np.diff(np.insert(indices, 0, -2)) > 1)[0]
    motif_loc = indices[start_indices]

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

    # sub1 = ts[loc1_st : loc1_end + m ]
    # for jj in range(loc2_st , loc2_end):
    #     sub2 = ts[jj : jj + m]
    #     temp_mp = mpx.compute(sub1,m, sub2)['mp']
    #     temp_bsf , temp_loc = bsfMotif(temp_mp)
    #     if temp_bsf <= bsf:
    #         bsf = temp_bsf
    #         bsf_loc = [temp_loc[0] + loc1_st, jj]

    indices = np.unique(np.concatenate((np.arange(loc1_st, loc1_end+m), np.arange(loc2_st, loc2_end+m))))
    temp_mp = mpx.compute(np.concatenate((ts[loc1_st : min(loc1_end + m, loc2_st) ], ts[loc2_st : loc2_end + m])),m)['mp']
    temp_bsf , temp_bsf_loc = bsfMotif(temp_mp)
    if temp_bsf <= bsf:
        bsf = temp_bsf
        bsf_loc = [indices[temp_bsf_loc[0]], indices[temp_bsf_loc[1]]]

    return bsf, bsf_loc, local_bsf


def prune(T, m, absf, bsf, amp, idxList):
    idx = findIndices(amp, m, absf, bsf, len(T))
    pruned_ts, pruned_idxList = T[idx], idxList[idx]
    return pruned_ts, pruned_idxList


def findIndices(amp, m, absf, bsf, n):
    # print('amp in pruning: ', amp)
    initset = np.where((amp >= absf) & (amp <= bsf))[0]
    # print(initset)
    split_indices = np.where(np.diff(initset) != 1)[0] + 1
    # print(split_indices)
    subarrays = np.split(initset, split_indices)
    margin = m //4

    for ii in range(len(subarrays)):
        lastitem = subarrays[ii][-1]
        if ii < len(subarrays) - 1:
            firstitem = subarrays[ii+1][0]
        else:
            firstitem = n
        gap = firstitem  - lastitem
        if gap < m+margin:
            subarrays[ii] = np.concatenate((subarrays[ii], np.arange(lastitem+1, lastitem + gap)))
        else:
            subarrays[ii] = np.concatenate((subarrays[ii], np.arange(lastitem+1, lastitem + m + margin)))


        firstitem = subarrays[ii][0]
        if ii > 0:
            lastitem = subarrays[ii-1][-1]
        else:
            lastitem = 0
        gap = firstitem  - lastitem
       
        if gap < margin:
            subarrays[ii] = np.concatenate((np.arange(firstitem - gap, firstitem), subarrays[ii]))
        else:
            subarrays[ii] = np.concatenate((np.arange(firstitem - margin, firstitem), subarrays[ii]))

    candidates = np.concatenate(subarrays)
    # candidates = np.unique(candidates)

    return candidates



if __name__=="__main__":
    T = genData(128, 8, dtype = 'rwalk')
    m = 8
    approxMP(T, m, dsr = 2)