function converted = gprof_convert(query, organism, target_namespace)
arguments
    query
    organism = 'rnorvegicus'
    target_namespace = 'ENSG'
end

query=cellstr(query(:)'); %is this enough to handle string vector?

args=pyargs('base_url','https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15/');
gp=py.gprofiler.GProfiler(args);
pyquery=py.list(cellstr(query)); %needs a row-vector cell array of chars
tic
gConvArgs=pyargs('organism',organism,'query',pyquery,'target_namespace',target_namespace);
res=gp.convert(gConvArgs); %returns list of dicts
res=Core_py2matlab(res); 
res=[res{:}];
converted=struct2table(res);
toc
