% Generative model
% 
% Author: David Meijer, 2021
% Refactoring and Docs: Roberto Barumerli, 2021
% Acoustic Research Institute, Wien
%
% DESC
% generate a sequence corrupted by noise sampled from the dist. p(x_t|s_t)
%
% IN 
% T         (int)   sequence total length in time units
% tchange   (int)   change point time index
% sigma     (float) standard deviation of random process with zero mean
% N         (int)   repetitions
%
% OUT
% x     (NxT float) N sequences with length T corrupted by noise
% s     (NxT float) N sequences with length T

function [x,s] = generativeModel(T,t_change,sigma,N)
    assert(t_change <= T, ...
        't_change is greater than the length of the sequence T')
    s = [-1*ones(N,t_change-1), ones(N,T-t_change+1)];
    x = s+sigma*randn(N,T);
end