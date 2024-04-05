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
        data = rnd.randn(n)/10
        data[30000:30000+len(motif)] = motif
        data[42000:42000+len(motif)] = motif
        data = zscore(data + rnd.randn(len(data))/20)
        

    if dtype == 'rwalk':
        # n = 64000
        # t = np.arange(0, 10, 0.01)
        # motif = np.abs(np.sin(t))
        # data = zscore(generate_random_walk(n))
        # data[30000:30000+len(motif)] = motif
        # data[42000:42000+len(motif)] = motif
        # data = zscore(data + rnd.randn(len(data))/20)
        data = generate_random_walk(n) /10
        data[10000:10300]  = np.abs(np.sin(np.arange(0,30,0.1))) /4
        data[40000:40300]  = np.abs(np.sin(np.arange(0,30,0.1)))/ 4
        data = zscore(data + np.random.rand(len(data))/20)

    if dtype == 'test':
        t = np.arange(0, 30, 1)
        motif = np.abs(np.sin(t))
        data = rnd.randn(n)/10
        data[20:20+len(motif)] = motif
        data[100:100+len(motif)] = motif
        data = zscore(data + rnd.randn(len(data))/20)
    
    return data


def generate_random_walk(n):
    np.random.seed(0)
    # # Generate random increments (+1 or -1) for each step
    # increments = np.random.choice([-1, 1], size=n-1)
    
    # # Calculate the cumulative sum of the increments to get the random walk
    # random_walk = np.cumsum(increments)
    
    # # Add starting point (0) to the random walk
    # random_walk = np.insert(random_walk, 0, 0)
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
