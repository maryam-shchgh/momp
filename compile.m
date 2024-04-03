rng(0);
clc;
n = 20000; m = 500;
T = randn(n,1)/10;
motif = sin(1:0.1:100)/20;

T(2000: 2000 + length(motif) -1) = motif;
T(6800: 6800 + length(motif) -1) = motif;
T = zscore(T + randn(length(T),1) /20);
length(motif)
[momp_out, mp_out] = momp(T, m, 1);

% [mp, ~, motifsIdx, ~] = mpx_v2(T, m, m);
% plot(mp)
% display(motifsIdx(1:2,1))