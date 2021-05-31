function [go_genes, go_geneTab] = get_go_genes(go_ids, organism, target_namespace, method)
arguments
    go_ids
    organism = 'rnorvegicus'
    target_namespace = 'ENSG'
    method='gprof'
end

switch method
    case 'gprof'
        gp=py.gprofiler.GProfiler();
        ids=py.list(cellstr(go_ids(:)')); %needs a row-vector cell array of chars
        gConvArgs=pyargs('organism',organism,'query',ids,'target_namespace',target_namespace);
        res=gp.convert(gConvArgs); %returns list of dicts
        res=Core_py2matlab(res); 
        res=[res{:}];
        go_geneTab=struct2table(res); %all IDs map, some gene symbols change.
        go_genes=unique(go_geneTab(:,{'converted','name'}),'rows');
        go_genes=renamevars(go_genes,'converted','id');
        go_genes=sortrows(go_genes,'name');
    otherwise
        error('unknown method')
end