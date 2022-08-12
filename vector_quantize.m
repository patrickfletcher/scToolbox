function result = vector_quantize(X, k, options, statopts) %, methodpars
arguments
    X
    k = floor(sqrt(size(X,1)))
    options.method = 'kmeans'
    options.Distance = 'sqeuclidean'
    options.Replicates = 1
    options.OnlinePhase='off'
    statopts.UseParallel = true
    statopts.MaxIter = 100
    statopts.Display = 'off'
end
% perform vector quantization of observations (rows) in X

statopts=statset(statopts);

switch options.method
    case 'kmeans'
    [idx,C, sumD] = kmeans(X, k, ...
        Distance=options.Distance, ...
        OnlinePhase=options.OnlinePhase, ...
        Replicates=options.Replicates, ...
        Options=statopts);

    case 'kmedoids'
    [idx,C, sumD] = kmedoids(X, k, ...
        Distance=options.Distance, ...
        OnlinePhase=options.OnlinePhase, ...
        Replicates=options.Replicates, ...
        Options=statopts);
end

result.k=k;
result.idx=idx;
result.C=C;
result.sumD=sumD;