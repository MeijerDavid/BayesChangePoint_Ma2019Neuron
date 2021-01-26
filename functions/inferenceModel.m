%Inference model
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