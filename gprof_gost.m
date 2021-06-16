function res = gprof_gost(id_list, organism, options)
arguments
    id_list
    organism
    options.Ordered (1,1) {mustBeNumericOrLogical} = false
    options.Sources = {'GO'}
    options.no_evidences = false
    options.user_threshold = 0.05
    options.significance_threshold_method = 'g_SCS'
end

% g:Profiler query
gp=py.gprofiler.GProfiler();

query=py.list(cellstr(id_list(:))'); %must be a row vector
sources=py.list(options.Sources);
gpargs=pyargs('organism',organism,'ordered',true,...
    'query',query,'sources',sources,'no_evidences',false,...
    'user_threshold',options.user_threshold,...
    'significance_threshold_method',options.significance_threshold_method);
res=gp.profile(gpargs);

res=Core_py2matlab(res);
res=[res{:}];
res=struct2table(res);