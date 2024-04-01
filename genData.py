import numpy as np
from numpy import random as rnd
from scipy.stats import zscore


def genData():
    np.random.seed(0)
    data = rnd.randn(64000)/10
    t = np.arange(0, 10, 0.01)
    motif = np.abs(np.sin(t))/10
    motif = motif + rnd.randn(len(motif))/20
    data[30000:30000+len(motif)] = motif
    data[42000:42000+len(motif)] = motif
    data = zscore(data + rnd.randn(len(data))/20)
    
    return data