% Generate responses
% 
% Author: David Meijer, 2021
% Refactoring and Docs: Roberto Barumerli, 2021
% Acoustic Research Institute, Wien
%
% DESC
% simulate estimation process N times
%
% IN 
% T         (int)   sequence total length in time units
% tchange   (int)   change point time index
% sigma     (float) standard deviation of random process with zero mean
% N         (int)   repetitions
%
% OUT
% responses (Nx1 int) index of the estimated t_change for N trials


function responses = genResponses(T,t_change,sigma,N)
    x = generativeModel(T,t_change,sigma,N);
    responses = inferenceModel(x,sigma);
end