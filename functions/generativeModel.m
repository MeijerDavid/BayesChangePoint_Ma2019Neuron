%Generative model
function [x,s] = generativeModel(T,tchange,sigma,N)
    s = [-1*ones(N,tchange-1), ones(N,T-tchange+1)];
    x = s+sigma*randn(N,T);
end