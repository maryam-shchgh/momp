import numpy as np
from numpy import random as rnd
from scipy.stats import zscore
import matplotlib.pyplot as plt

# %matplotlib widget

def genData(n, dtype = 'rwalk'):
    np.random.seed(0)
    
    if dtype == 'random':
        n = 64000
        t = np.arange(0, 10, 0.01)
        motif = np.abs(np.sin(t))
        T = rnd.randn(n)/10
        T[30000:30000+len(motif)] = motif
        T[42000:42000+len(motif)] = motif
        T = zscore(T + rnd.randn(len(T))/20)
        

    if dtype == 'rwalk':
        T = generate_random_walk(n) /10
        T[30000:30600]  = np.abs(np.sin(np.arange(0,6,0.01))) /4
        T[60000:60600]  = np.abs(np.sin(np.arange(0,6,0.01)))/ 4
        T = zscore(T + np.random.rand(len(T))/20)

    if dtype == 'test':
        t = np.arange(0, 30, 1)
        motif = np.abs(np.sin(t))
        T = rnd.randn(n)/10
        T[20:20+len(motif)] = motif
        T[100:100+len(motif)] = motif
        T = zscore(T + rnd.randn(len(T))/20)
    
    return T


def generate_random_walk(n):
    np.random.seed(0)
    # Probability to move up or down
    prob = [0.05, 0.95]  
    
    # statically defining the starting position
    start = 2 
    random_walk = [start]
    
    # creating the random points
    rr = np.random.random(n)
    downp = rr < prob[0]
    upp = rr > prob[1]
    
    
    for idownp, iupp in zip(downp, upp):
        down = idownp and random_walk[-1] > 1
        up = iupp and random_walk[-1] < 4
        random_walk.append(random_walk[-1] - down + up)
    
    random_walk = np.array(random_walk)
    return random_walk
