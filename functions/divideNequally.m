%Divide N equally over n places (some randomly get at most 1 more)
function counts = divideNequally(N,n)
    base = floor(N/n);
    counts = base*ones(1,n);
    idx_1extra = randperm(n,N-base*n);
    counts(idx_1extra) = counts(idx_1extra)+1;
end