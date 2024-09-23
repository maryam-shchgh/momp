# iMOMP : interleaved MOMP 

**Use iMOMP for Motif Discovery in a Set of Timeseries**

## Table of Contents
1. [Installation](#installation)
2. [How to run iMOMP?](#usage)

## Installation

1. **Clone the repository:**
   git clone https://github.com/maryam-shchgh/momp.git

2. **Run Usage:**
   Following the provided usage, you can easily study your set of timeseries.
   Current version of iMOMP is in MATLAB.

## usage
1. [momp_out, momp_loc] = mist_v02(Ts_set, m, verbose, run_mpx);
2. **Parameters:**
   Ts_set: Timeseries set, with each Ts organized in a separate column.
   m: Motif length.
   verbose (optional): Set to 1 to enable terminal logging (default = 1).
   run_mpx (optional): Set to 1 to run Mpx on all timeseries (default = 0).
