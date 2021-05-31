function res = gProfileQuery(id_list, organism)
arguments
    id_list
    organism
    options.Ordered (1,1) {mustBeNumericOrLogical} = false
    options.Sources {mustBe
    options.
end

% g:Profiler query
gp=py.gprofiler.GProfiler();

query=py.list(cellstr(id_list(:))'); %must be a row vector
sources=py.list({'GO:MF','GO:CC','GO:BP'});
gpargs=pyargs('organism','rnorvegicus','ordered',true,...
    'query',query,'sources',sources);
res=gp.profile(gpargs);

res=Core_py2matlab(res);
res=[res{:}];