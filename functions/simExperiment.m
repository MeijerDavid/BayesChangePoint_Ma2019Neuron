%Simulate an experiment for 1 subject
function [responses,true_tchanges] = simExperiment(T,sigma,N)
    
    %Assume that tchange can take any position between 1 and T
    %Divide the number of trials per tchange equally.
    nTrials_tchange = divideNequally(N,T);
    
    true_tchanges = nan(N,1);
    responses = nan(N,1);
    idx_begin = 1;
    for tchange=1:T
        idx_end = idx_begin+nTrials_tchange(tchange)-1;
        true_tchanges(idx_begin:idx_end) = tchange;
        responses(idx_begin:idx_end) = genResponses(T,tchange,sigma,nTrials_tchange(tchange));
        idx_begin = idx_begin+nTrials_tchange(tchange);
    end
    
    %Randomize the trial order (as in real experiment)
    random_order = randperm(N);
    responses = responses(random_order);
    true_tchanges = true_tchanges(random_order);
end