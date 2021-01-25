%Compute the log-likelihood
function LL = compLogLike(log_sigma,responses,true_tchanges,T,N)
    
    %Use log_sigma because it is unbounded, whereas sigma has to be > 0
    sigma = exp(log_sigma);

    all_tchange = unique(true_tchanges);
    n_tchange = length(all_tchange);
    
    bin_edges = [0.5, (1:T)+0.5];
    
    loglikes = nan(1,n_tchange);
    for i=1:n_tchange
        
        resp_subj = responses(true_tchanges==all_tchange(i));
        counts_subj = histcounts(resp_subj,bin_edges);
        
        resp_sim = genResponses(T,all_tchange(i),sigma,N);
        predicted_prob = histcounts(resp_sim,bin_edges)/N + eps;  
        %+eps is hack to avoid issues: log(0)=-inf and 0*log(0)=NaN
        
        loglikes(i) = sum(counts_subj.*log(predicted_prob));             
    end
    LL = sum(loglikes);
end