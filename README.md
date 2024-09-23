# iMOMP : interleaved MOMP 

**Use iMOMP for Motif Discovery in a Set of Timeseries**

## Table of Contents
1. [Installation](#installation)
2. [How to run iMOMP?](#usage)
3. [Related Paper(s)](#references)

## Installation

1. **Clone the repository:**
   git clone https://github.com/maryam-shchgh/momp.git

2. **Run Usage:**
   Following the provided usage, you can easily study your set of timeseries.
   Current version of iMOMP is in MATLAB.

## How to run iMOMP
1. **Usage:**
   [momp_out, momp_loc] = mist_v02(Ts_set, m, verbose, run_mpx);
3. **Parameters:**
   Ts_set: Timeseries set, with each Ts organized in a separate column.
   m: Motif length.
   verbose (optional): Set to 1 to enable terminal logging (default = 1).
   run_mpx (optional): Set to 1 to run Mpx on all timeseries (default = 0).

## Related Paper(s)
This code is associated with the paper titled:

**Title:** Matrix Profile XXXI: Motif-Only Matrix Profile: Orders of Magnitude Faster  
**Authors:** Maryam Shahcheraghi, Chin Chia Michael Yeh, Yan Zheng, Junpeng Wang, Zhongfang Zhuang, Liang Wang, Mahashweta Das, and Eamonn Keogh
**Published in:** IEEE ICKG 2024  
**Link:** [Read the full version of the Paper](https://lnkd.in/guwmUJaf)

### Abstract
Approximately repeated subsequences in a longer timeseries, i.e.,time series motifs, are important primitive in time series data mining. Motifs are used in dozens of downstream tasks, including classification, clustering, summarization, rule discovery, segmentation etc. Time series motif discovery is a notoriously computationally expensive task. Some motif discovery algorithms are fast in the best case, but in other datasets, even if both the data and motif lengths are held the same, both their time and space complexity can explode. The Matrix Profile has the nice property that its time and space complexity are independent of the data. Moreover, the Matrix Profile is fast enough for datasets with about one million datapoints, which covers a large fraction of user cases. However, there are situations where we may wish to consider datasets which are much larger. In this work, we introduce the first lower bound for the Matrix Profile and an algorithm that exploits that lower bound to allow orders of magnitude speed up for exact motif search on real-world datasets. We demonstrate the utility of our ideas with the largest and most ambitious motif discovery experiments ever attempted.

### Key Contributions
- Introduction of the lbMP lower bound for the Matrix Profile.
- Development of the MOMP algorithm for exact motif discovery by up to two orders of magnitude faster.

### Citation
If you use this work in your research, please cite it as follows: [The link will be inserted soon]
