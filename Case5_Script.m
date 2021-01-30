% Play with Bayesian Change-Point Detection Model by Wei Ji Ma (2019) 
% https://doi.org/10.1016/j.neuron.2019.09.037
% 
% Author: David Meijer, 2021
% Refactoring and Docs: Roberto Barumerli, 2021
% Acoustic Research Institute, Wien

clear var
close all

% Add functions path
addpath(genpath('functions'));

% number of repetitions - montecarlo simulations
num_exp = 10e3;

%% Question 3.5.i - generative model and inference model
sequence_length = 10; 
sigma = 1; 

prop_correct = nan(1,sequence_length);
for t_change=1:sequence_length
    responses = genResponses(sequence_length,t_change,sigma,num_exp);
    prop_correct(t_change) = sum(responses==t_change)/num_exp;
end
fprintf('Computed mean of correct answers: %0.2f\n', ...
                mean(prop_correct))

figure; 
plot(1:sequence_length,prop_correct,'k-o');
ylim([0.5 1]); 
xlabel('True time of change [#]'); 
ylabel('Proportion correct');
title('3.5.i - ratio of corrected answers with different t-changes');

%% Question 3.5.k - Exploration of model parameters (sigma and T)
sequence_length = 2:2:16;
sigma = 1:3; 

mean_prop_correct = zeros(length(sigma),length(sequence_length));
responses = zeros(num_exp,1);

for i = 1:length(sigma)
    for j=1:length(sequence_length)
        for t_change=1:sequence_length(j)
            responses = genResponses(sequence_length(j),t_change,sigma(i),num_exp);
            mean_prop_correct(i,j) = mean_prop_correct(i,j) + mean(sum(responses==t_change)/num_exp);
        end
        mean_prop_correct(i,j) = mean_prop_correct(i,j)/sequence_length(j);
    end
end

figure; 
hold on; 
colors = [[0 0 1]; [1 0 0]; [.8 .8 0]]; 

for i=1:length(sigma)
    plot(sequence_length,mean_prop_correct(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ', num2str(sigma(i)), ' deg'];
end
ylim([0.2 1]); xlabel('Sequence length (T)'); ylabel('Proportion correct');    
legend(leg_labels,'location','northeast','Interpreter','latex');
title('3.5.k - Proportion correct for different sequence lengths'); 

%% Parameter recovery (this takes about 1 minute to run)
%As a function of the number of behavioral trials
sequence_length = 10;
sigma = 1:3; 
num_trials = [10, 20, 50, 100, 200, 500, 1000];
num_repeats = 100; Nsims = 1000;
fitted_sigma = nan(length(sigma),length(num_trials),num_repeats);
for i = 1:(length(sigma)*length(num_trials))
    disp(['Starting ' num2str(i) ' of ' num2str(length(sigma)*length(num_trials)) ' with ' num2str(num_repeats) ' each.']); 
    [i_sigma,j_nTrials] = ind2sub([length(sigma),length(num_trials)],i);
    for i_repeat=1:num_repeats
        %Simulate an experiment for 1 subject
        [responses,true_tchanges] = simExperiment(sequence_length,sigma(i_sigma),num_trials(j_nTrials));
        
        %Negative log likelihood is the objective function for fitting
        fitfun = @(log_sigma) -compLogLike(log_sigma,responses,true_tchanges,sequence_length,Nsims);
        
        %Set min-max bounds for the (log-)sigma as an aid to the fitting algorithm
        log_sigma_bounds = [log(0.1),log(20)];
        
        %Fit the model to this subject's dataset and obtain a fitted sigma
        fitted_sigma(i_sigma,j_nTrials,i_repeat) = exp(fminbnd(fitfun,log_sigma_bounds(1),log_sigma_bounds(2)));
    end
end
mean_fitted_sigma = mean(fitted_sigma,3);
std_fitted_sigma = std(fitted_sigma,[],3);

figure; colors = [[0 0 1]; [1 0 0]; [.8 .8 0]]; h = nan(1,3); 
subplot(1,2,1); hold on; title('Fitted Sigma Means');
for i=1:length(sigma)
    h(i) = plot(num_trials,mean_fitted_sigma(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ' num2str(sigma(i)) ' deg'];
end
ylim([0 4]); xlabel('nTrials'); ylabel('Mean of fitted sigma ( deg)');    
legend(h,leg_labels,'location','northeast');

subplot(1,2,2); hold on; title('Fitted Sigma SDs');
for i=1:length(sigma)
    h(i) = plot(num_trials,std_fitted_sigma(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ' num2str(sigma(i)) ' deg'];
end
ylim([0 3]); xlabel('nTrials'); ylabel('SD of fitted sigma ( deg)');    
legend(h,leg_labels,'location','northeast');
%It is clear that you get decent estimates of sigma with at least 100
%trials total (i.e. 10 trials per t_change). Fewer trials lead to biased
%estimates, whereas more trials only marginally reduce the variability of 
%the fitted sigmas.

%% Power analysis (this takes several hours to run)
%As a function of effect size and group size.
sequence_length = 10; num_trials = 100; num_repeats = 200; Nsims = 1000;
effect_sizes = 0.1:0.2:0.9;                                                 %Cohen's d.
num_subjs = 10:10:50;
log_sigma_mean = log(2); log_sigma_SD = log(0.5);                           %For sampling behavioral sigmas
log_sigma_bounds = [log(0.1),log(20)];                                      %For the fitting algorithm

prop_sign = nan(length(num_subjs),length(effect_sizes));
for i = 1:(length(num_subjs)*length(effect_sizes))
    [i_numsubj,j_effsize] = ind2sub([length(num_subjs),length(effect_sizes)],i);
    disp(['Starting ' num2str(i) ' of ' num2str(length(num_subjs)*length(effect_sizes)) ' with effect size = ' num2str(effect_sizes(j_effsize)) ' and ' num2str(num_subjs(i_numsubj)) ' subjects.']); 
    
    %Translate effect size into difference of log_sigmas
    log_sigma_diff = effect_sizes(j_effsize)*log_sigma_SD;
    
    H = nan(1,num_repeats); %initialize
    
    %Repeat a number of simulated experiments (each experiment has multiple subjects)   
    for i_repeat=1:num_repeats
        
        %Sample behavioral sigmas from a normal distribution (entire group)   
        real_log_sigmas_1 = log_sigma_mean + log_sigma_SD*randn(1,num_subjs(i_numsubj));
        real_log_sigmas_2 = real_log_sigmas_1 + log_sigma_diff;                         
        
        fitted_log_sigmas_1 = nan(size(real_log_sigmas_1));
        fitted_log_sigmas_2 = nan(size(real_log_sigmas_2));
        
        %Simulate an experiment per subject and fit the data
        for i_subj = 1:num_subjs(i_numsubj)
            
            %Condition 1
            [responses,true_tchanges] = simExperiment(sequence_length,real_log_sigmas_1(i_subj),num_trials);
            fitfun = @(log_sigma) -compLogLike(log_sigma,responses,true_tchanges,sequence_length,Nsims);
            fitted_log_sigmas_1(i_subj) = fminbnd(fitfun,log_sigma_bounds(1),log_sigma_bounds(2));
            
            %Condition 2
            [responses,true_tchanges] = simExperiment(sequence_length,real_log_sigmas_2(i_subj),num_trials);
            fitfun = @(log_sigma) -compLogLike(log_sigma,responses,true_tchanges,sequence_length,Nsims);
            fitted_log_sigmas_2(i_subj) = fminbnd(fitfun,log_sigma_bounds(1),log_sigma_bounds(2));
        end
        
        %Perform paired t-test to see if it is significant (at alpha = 0.05)   
        H(i_repeat) = ttest(fitted_log_sigmas_1,fitted_log_sigmas_2);
    end
    
    %Proportion of experiments that reached statistical significance
    prop_sign(i_numsubj,j_effsize) = sum(H)/num_repeats;
end

figure;
imagesc(effect_sizes,num_subjs,prop_sign,[0 1]); 
set(gca,'YDir','normal'); title('Power analysis');
set(gca,'XTick',effect_sizes,'XTickLabel',effect_sizes); xlabel('Effect size (Cohen''s d)');
set(gca,'YTick',num_subjs,'YTickLabel',num_subjs); ylabel('Number of subjects');
h = colorbar; ylabel(h,'power (proportion significant ttests');
hold on; contour(effect_sizes,num_subjs,prop_sign,[.8 .8],'LineColor','r','Linewidth',2,'ShowText','on'); 

%Somewhat surprisingly, power decreases with effect sizes > 0.5. 
%I believe that this happens because with larger effect size the sensory
%noise of the 2nd condition also gets larger. At such large sensory noise 
%levels the other experimental settings (e.g. spatial locations -1 and +1,
%100 behavioral trials only) may not allow for accurate estimates of these
%large sigmas, thus leading to additional noise in the data..
