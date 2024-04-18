rng(0);
clc;
verbose = 1;
n = 2000; m = 64;

T = genData(n);

momp(T, m, verbose);
