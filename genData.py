import numpy as np
from numpy import random as rnd
from scipy.stats import zscore
import matplotlib.pyplot as plt

# %matplotlib widget

def genData(dtype = 'rwalk'):
    n = 64000
    np.random.seed(0)
    t = np.arange(0, 10, 0.01)
    motif = np.abs(np.sin(t))
    # motif = motif + rnd.randn(len(motif))/20

    if dtype == 'random':
        data = rnd.randn(n)/10
        # t = np.arange(0, 10, 0.01)
        # motif = np.abs(np.sin(t))/10
        # motif = motif + rnd.randn(len(motif))/20
        

    if dtype == 'rwalk':
        data = zscore(generate_random_walk(n))

    data[30000:30000+len(motif)] = motif
    data[42000:42000+len(motif)] = motif
    data = zscore(data + rnd.randn(len(data))/20)
    
    return data


def generate_random_walk(n):
    # Generate random increments (+1 or -1) for each step
    increments = np.random.choice([-1, 1], size=n-1)
    
    # Calculate the cumulative sum of the increments to get the random walk
    random_walk = np.cumsum(increments)
    
    # Add starting point (0) to the random walk
    random_walk = np.insert(random_walk, 0, 0)
    
    return random_walk




if __name__=="__main__":
    genData(n=64000, dtype = 'rwalk')