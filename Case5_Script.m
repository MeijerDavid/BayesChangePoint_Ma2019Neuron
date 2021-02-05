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

%% Question 3.5.i - generative model
% sequence_length = 10; 
% sigma = 1; 
% 
% prop_correct = nan(1,sequence_length);
% for t_change=1:sequence_length
%     responses = genResponses(sequence_length,t_change,sigma,num_exp);
%     prop_correct(t_change) = sum(responses==t_change)/num_exp;
% end
% fprintf('Computed mean of correct answers: %0.2f\n', ...
%                 mean(prop_correct))
% 
% figure; 
% plot(1:sequence_length,prop_correct,'k-o');
% ylim([0.5 1]); 
% xlabel('True time of change [#]'); 
% ylabel('Proportion correct');
% title('3.5.i - ratio of corrected answers with different t_changes');

%% Question 3.5.k - sensibility analysis of the parameters
% sequence_length = 2:2:16;
% sigma = 1:3; 
% 
% mean_prop_correct = zeros(length(sigma),length(sequence_length));
% responses = zeros(N,1);
% 
% for i = 1:length(sigma)
%     for j=1:length(sequence_length)
%         for t_change=1:sequence_length(j)
%             responses = genResponses(sequence_length(j),t_change,sigma(i),num_exp);
%             mean_prop_correct(i,j) = mean_prop_correct(i,j) + ...
%                 mean(sum(responses==t_change)/num_exp);
%         end
%         mean_prop_correct(i,j) = mean_prop_correct(i,j)/sequence_length(j);
%     end
% end
% 
% figure; 
% hold on; 
% colors = [[0 0 1]; [1 0 0]; [.8 .8 0]]; 
% 
% for i=1:length(sigma)
%     plot(sequence_length,mean_prop_correct(i,:),'-o','Color',colors(i,:)); 
%     leg_labels{i} = ['sigma = ', num2str(sigma(i)), ' deg'];
% end
% ylim([0.2 1]); xlabel('T'); ylabel('Proportion correct');    
% legend(leg_labels,'location','northeast','Interpreter','latex');
% title('3.5.k - ratio of corrected answers with different sequence lengths'); 

%% Parameter recovery (this takes about 20 minutes to run)
sequence_length = 10;
sigma = 1:3; 
num_trials = [10 20, 50, 100, 200, 500, 1000];
num_repeats = 100; 
num_sims = 1000; % montecarlo on parameter estimation

fitted_sigma = zeros(length(sigma),length(num_trials),num_repeats);

for i = 1:(length(sigma)*length(num_trials))
    disp(['Starting ' num2str(i) ' of ' num2str(length(sigma)*length(num_trials)) ' with ' num2str(num_repeats) ' repetitions each.']); 
    [i_sigma,j_nTrials] = ind2sub([length(sigma),length(num_trials)],i);
    for i_repeat=1:num_repeats
        % Simulate an experiment for 1 subject
        [responses, t_changes] = ...
            simExperiment(sequence_length, ...
                            sigma(i_sigma), ...
                            num_trials(j_nTrials));
        
        % Negative log likelihood is the objective function for fitting
        fitfun = @(log_sigma) -compLogLike(log_sigma, ...
                                            responses, ...
                                            t_changes, ...
                                            sequence_length, ...
                                            num_sims);
        
        % Randomly select a starting value for free parameter: log(sigma)
        % Select from uniform distribution between log(0.5) and log(3.5)
        log_sigma_start = rand(1)*(log(3.5)-log(0.5))+log(0.5);
        
        % Fit the model to this subject's dataset and obtain a fitted sigma
        fitted_sigma(i_sigma,j_nTrials,i_repeat) = exp(fminsearch(fitfun,log_sigma_start));
        % Many attempts did not converge (more evaluations needed). But this will do for now.
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
ylim([0 3]); xlabel('nTrials'); ylabel('Mean of fitted sigma ( deg)');    
legend(h,leg_labels,'location','northeast');

subplot(1,2,2); hold on; title('Fitted Sigma SDs');
for i=1:length(sigma)
    h(i) = plot(num_trials,std_fitted_sigma(i,:),'-o','Color',colors(i,:)); 
    leg_labels{i} = ['sigma = ' num2str(sigma(i)) ' deg'];
end
ylim([0 3]); xlabel('nTrials'); ylabel('SD of fitted sigma ( deg)');    
legend(h,leg_labels,'location','northeast');
%I imagined these results looking a little nicer. Might have to do with the
%bad methods here. E.g. not using multiple starting points, not having
%fminsearch converge, etc.

%% Power analysis
% Not done because previous results were already so messy ..
% and playing time's up! Time for weekend :)
