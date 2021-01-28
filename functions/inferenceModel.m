% Inference model
% 
% Author: David Meijer, 2021
% Refactoring and Docs: Roberto Barumerli, 2021
% Acoustic Research Institute, Wien
%
% DESC
% Return an optimal estimate of t_change
%
% IN 
% x         (NxT float) N sequences with length T corrupted by noise
% sigma     (float) standard deviation of random process with zero mean
%
% OUT
% response  (int) estimated index of t_change
% MAP       (float) posterior probability of the answer

function [response,MAP] = inferenceModel(x,sigma)
    if nargout == 2
        likelihood = exp((2/sigma^2)*cumsum(x,2,'reverse'));
        posterior = bsxfun(@rdivide,likelihood,sum(likelihood,2));
        [MAP,response] = max(posterior,[],2);
    else
        %The above is correct. But this is faster.
        [~,response] = max(cumsum(x,2,'reverse'),[],2);
        MAP = 'dummy';
    end
    
end