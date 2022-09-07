function converted = gprof_convert(query, organism, target_namespace, options)
arguments
    query
    organism = 'rnorvegicus'
    target_namespace = 'ENSG'
    options.base_url = []
end

query=cellstr(query(:)'); %is this enough to handle string vector?

% if ~isempty(options.base_url)
args=pyargs('base_url','https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15');
gp=py.gprofiler.GProfiler(args);
% gp=py.gprofiler.GProfiler();

pyquery=py.list(cellstr(query)); %needs a row-vector cell array of chars
tic
gConvArgs=pyargs('organism',organism,'query',pyquery,'target_namespace',target_namespace);
res=gp.convert(gConvArgs); %returns list of dicts
res=Core_py2matlab(res); 
res=[res{:}];
converted=struct2table(res);
toc
