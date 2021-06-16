function converted = gprof_convert(query, organism, target_namespace)
arguments
    query
    organism = 'rnorvegicus'
    target_namespace = 'ENSG'
end

query=cellstr(query(:)'); %is this enough to handle string vector?

gp=py.gprofiler.GProfiler();
pyquery=py.list(cellstr(query)); %needs a row-vector cell array of chars
tic
gConvArgs=pyargs('organism',organism,'query',pyquery,'target_namespace',target_namespace);
res=gp.convert(gConvArgs); %returns list of dicts
res=Core_py2matlab(res); 
res=[res{:}];
converted=struct2table(res);
toc
