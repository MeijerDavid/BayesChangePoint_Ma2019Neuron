% Author: David Meijer, 2021
% Refactoring and Docs: Roberto Barumerli, 2021
% Acoustic Research Institute, Wien
%
% DESC
% Compute the log-likelihood of experimental data, is
%
% IN 
% responses (Nx1 int) index of the estimated t_change for N trials
% t_changes (Nx1 int) index of the real t_change for N trials
% T         (int)   sequence total length in time units
% N         (int)   repetitions
%
% OUT
% LL (float) log-likelihood of the value log_sigma


function LL = compLogLike(log_sigma,responses,t_changes,T,N)
    %Use log_sigma because it is unbounded, whereas sigma has to be > 0
    sigma = exp(log_sigma);
    
    assert(T == length(unique(t_changes)), "something wrong with the data provided")
    
    bin_edges = [0.5, (1:T)+0.5];
    
    loglikes = zeros(1,T);
    
    for t=1:T % cardinality of t_changes (i is also the value of t_change)
        resp_subj = responses(t_changes == t); % select trials with specific t_change
        counts_subj = histcounts(resp_subj, bin_edges); % real-data distribution
        
        resp_sim = genResponses(T, t, sigma, N); % simulate experiments with guessed sigma
        predicted_prob = histcounts(resp_sim,bin_edges)/N + eps; % simulated-data distribution
        %+eps is hack to avoid issues: log(0)=-inf and 0*log(0)=NaN
        
        loglikes(t) = sum(counts_subj.*log(predicted_prob));            
    end
    LL = sum(loglikes);
end