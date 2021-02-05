% Author: David Meijer, 2021
% Refactoring and Docs: Roberto Barumerli, 2021
% Acoustic Research Institute, Wien
%
% DESC
% Simulate an experiment for 1 subject
%
% IN 
% T         (int)   sequence total length in time units
% sigma     (float) standard deviation of random process with zero mean
% N         (int)   repetitions
%
% OUT
% responses (Nx1 int) index of the estimated t_change for N trials
% t_changes (Nx1 int) index of the real t_change for N trials


function [responses, t_changes] = simExperiment(T, sigma, N)
    % Assume that t_change can take any position between 1 and T
    % Divide the number of trials per t_change equally.
    nTrials_tchange = divideNequally(N,T);
    
    t_changes = zeros(N,1);
    responses = zeros(N,1);
    idx_begin = 1;
    for t=1:T
        idx_end = idx_begin+nTrials_tchange(t)-1;
        t_changes(idx_begin:idx_end) = t;
        responses(idx_begin:idx_end) = genResponses(T,t,sigma,nTrials_tchange(t));
        idx_begin = idx_begin+nTrials_tchange(t);
    end
    
    %Randomize the trial order (as in real experiment)
    random_order = randperm(N);
    responses = responses(random_order);
    t_changes = t_changes(random_order);
end