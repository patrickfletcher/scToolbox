function res = gprof_gost(id_list, organism, options)
arguments
    id_list
    organism
    options.Ordered = false
    options.Sources = {'GO'}
    options.no_evidences = false
    options.user_threshold = 0.05
    options.significance_threshold_method = 'g_SCS'
end

% g:Profiler query
args=pyargs('base_url','https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15/');
gp=py.gprofiler.GProfiler(args);

if isstring(id_list)||iscellstr(id_list)
    query=py.list(cellstr(id_list(:))'); %must be a row vector
else
    query=id_list; %assume it is a py.dict for multiquery - Ideally would build it here 
end

sources=py.list(options.Sources);
gpargs=pyargs('organism',organism,'ordered',options.Ordered,...
    'query',query,'sources',sources,'no_evidences',options.no_evidences,...
    'user_threshold',options.user_threshold,...
    'significance_threshold_method',options.significance_threshold_method);
res=gp.profile(gpargs);

res=Core_py2matlab(res);
res=[res{:}];
res=struct2table(res);

res.expectedfrac=res.term_size./res.effective_domain_size;
res.expected=res.expectedfrac.*res.query_size;
res.fold_enrichment=res.intersection_size./res.expected;