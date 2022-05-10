clear;
close all;
alpha=[-17.4, 12.7];   % True DOAs
SNR=10;
K=50;                  % The number of snapshots
M=10;                  % The number of sensors

resolution=6;          % grid interval
search_area=[-90:resolution:90];   % grid 
etc=M;                             %   etc may be chosen from 1 to M.

X=signal(M, alpha, SNR, K);     
[Pm_root,search_root]=Bayesian_DOA_root(X,search_area,etc);   %main
plot(search_root,Pm_root)
  