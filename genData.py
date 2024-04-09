import numpy as np
from numpy import random as rnd
from scipy.stats import zscore
import matplotlib.pyplot as plt

# %matplotlib widget

def genData(n, mlen, dtype = 'rwalk'):
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
        # ii, jj = round(0.3*len(T)) , round(0.7*len(T))
        ii, jj = 32, 100
        # T[ii : ii + mlen]  = np.abs(np.sin(np.arange(0,mlen*0.01,0.01))) /4
        T[ii : ii + mlen] = np.concatenate((np.arange(0,2,0.25), np.arange(2, 0, -0.25)))
        T[jj : jj + mlen]  = T[ii : ii + mlen]
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
    
    for _ in range(n-1):  # Loop n-1 times to generate a time series of length n
        rr = np.random.random()  # Generate a single random number
        if rr < prob[0] and random_walk[-1] > 1:
            random_walk.append(random_walk[-1] - 1)  # Move down
        elif rr > prob[1] and random_walk[-1] < 4:
            random_walk.append(random_walk[-1] + 1)  # Move up
        else:
            random_walk.append(random_walk[-1])  # Stay in the same position
    
    random_walk = np.array(random_walk)
    return random_walk
