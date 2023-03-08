function converted = gprof_convert(query, organism, target_namespace, options)
arguments
    query
    organism = 'rnorvegicus'
    target_namespace = 'ENSG'
    options.base_url = []
end

% if ~isempty(options.base_url)
args=pyargs('base_url','https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15');
gp=py.gprofiler.GProfiler(args);
% gp=py.gprofiler.GProfiler();

query=cellstr(query(:)'); %needs a row-vector cell array of chars
pyquery=py.list(query); 
gConvArgs=pyargs('organism',organism,'query',pyquery,'target_namespace',target_namespace);

tic
res=gp.convert(gConvArgs); %returns list of dicts
gpConvertTime = toc

tic
%convert the list of dicts to dict of lists in python:
res = pyrun("out={key: [i[key] for i in data] for key in data[0]}","out",data=res);
resStruct = struct(res);
fields=fieldnames(resStruct);
for i=1:length(fields)
    field = fields{i};
    thisData=resStruct.(field);
    pyType = class(thisData{1});
    switch pyType
        case {'py.str'}
            resStruct.(field) = string(thisData)';
        case {'py.int'}
            resStruct.(field) = double(thisData)';
    end
end
converted=struct2table(resStruct);
Python2MatlabTime = toc
end