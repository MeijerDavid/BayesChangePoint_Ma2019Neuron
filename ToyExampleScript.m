
%Play with Bayesian Change-Point Detection Model by Wei Ji Ma (2019)
%https://doi.org/10.1016/j.neuron.2019.09.037
clear all;

%Add functions path
addpath([cd filesep 'ToyExampleFunctions']);

%% Create figure 5.i
T = 10; sigma = 1; N = 10000;
prop_correct = nan(1,T);
for tchange=1:T
    responses = genResponses(T,tchange,sigma,N);
    prop_correct(tchange) = sum(responses==tchange)/N;
end
mean(prop_correct) %--> command window

figure; plot(1:T,prop_correct,'k-o'); ylim([0.5 1]); 
xlabel('True time of change'); ylabel('Proportion correct');
title('Figure 5.i (Ma 2019, Neuron)');

%% Create figure 5.k
T = 2:2:16; sigma = 1:3; N = 10000;
mean_prop_correct = nan(length(sigma),length(T));
for ind = 1:(length(sigma)*length(T))
    [i_sigma,j_T] = ind2sub([length(sigma),length(T)],ind);
    for tchange=1:T(j_T)
        responses = genResponses(T(j_T),tchange,sigma(i_sigma),N);
        mean_prop_correct(i_sigma,j_T) = mean(sum(responses==tchange)/N);
    end
end

figure; hold on; colors = [[0 0 1]; [1 0 0]; [.8 .8 0]]; h = nan(1,3);
for i=1:length(sigma)
    h(i) = plot(T,mean_prop_correct(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ' num2str(sigma(i)) '°'];
end
ylim([0.2 1]); xlabel('T'); ylabel('Proportion correct');    
legend(h,leg_labels,'location','northeast');
title('Figure 5.k (Ma 2019, Neuron)'); 
%I don't know why the figure looks different for sigma=2° and sigma=3°

%% Parameter recovery (this takes about 20 minutes to run)
T = 10; sigma = 1:3; nTrials = [10 20, 50, 100, 200, 500, 1000];
nRepeats = 100; Nsims = 1000;
fitted_sigma = nan(length(sigma),length(nTrials),nRepeats);
for ind = 1:(length(sigma)*length(nTrials))
    disp(['Starting ' num2str(ind) ' of ' num2str(length(sigma)*length(nTrials)) ' with ' num2str(nRepeats) ' each.']); 
    [i_sigma,j_nTrials] = ind2sub([length(sigma),length(nTrials)],ind);
    for i_repeat=1:nRepeats
        %Simulate an experiment for 1 subject
        [responses,true_tchanges] = simExperiment(T,sigma(i_sigma),nTrials(j_nTrials));
        
        %Negative log likelihood is the objective function for fitting
        fitfun = @(log_sigma) -compLogLike(log_sigma,responses,true_tchanges,T,Nsims);
        
        %Randomly select a starting value for free parameter: log(sigma)
        %Select from uniform distribution between log(0.5) and log(3.5)
        log_sigma_start = rand(1)*(log(3.5)-log(0.5))+log(0.5);
        
        %Fit the model to this subject's dataset and obtain a fitted sigma
        fitted_sigma(i_sigma,j_nTrials,i_repeat) = exp(fminsearch(fitfun,log_sigma_start));
        %Many attempts did not converge (more evaluations needed). But this will do for now.
    end
end
mean_fitted_sigma = mean(fitted_sigma,3);
std_fitted_sigma = std(fitted_sigma,[],3);

figure; colors = [[0 0 1]; [1 0 0]; [.8 .8 0]]; h = nan(1,3); 
subplot(1,2,1); hold on; title('Fitted Sigma Means');
for i=1:length(sigma)
    h(i) = plot(nTrials,mean_fitted_sigma(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ' num2str(sigma(i)) '°'];
end
ylim([0 3]); xlabel('nTrials'); ylabel('Mean of fitted sigma (°)');    
legend(h,leg_labels,'location','northeast');

subplot(1,2,2); hold on; title('Fitted Sigma SDs');
for i=1:length(sigma)
    h(i) = plot(nTrials,std_fitted_sigma(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ' num2str(sigma(i)) '°'];
end
ylim([0 3]); xlabel('nTrials'); ylabel('SD of fitted sigma (°)');    
legend(h,leg_labels,'location','northeast');
%I imagined these results looking a little nicer. Might have to do with the
%bad methods here. E.g. not using multiple starting points, not having
%fminsearch converge, etc.

%% Power analysis
%Not done because previous results were already so messy ..
%and playing time's up! Time for weekend :)
