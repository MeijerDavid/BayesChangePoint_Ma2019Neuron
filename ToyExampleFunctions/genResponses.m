%Generate responses
function responses = genResponses(T,tchange,sigma,N)
    x = generativeModel(T,tchange,sigma,N);
    responses = inferenceModel(x,sigma);
end