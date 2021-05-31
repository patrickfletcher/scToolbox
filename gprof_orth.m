function [orthgenes, orthTab] = gprof_orth(query, source, target)

gp=py.gprofiler.GProfiler();

pyquery=py.list(cellstr(query(:)'));
gOrthArgs=pyargs('organism',source,'query',pyquery,'target',target);

res=gp.orth(gOrthArgs); %returns list of dicts
res=Core_py2matlab(res); 
res=[res{:}];
orthTab=struct2table(res);

orthgenes=unique(orthTab(:,{'ortholog_ensg','name'}),'rows');
orthgenes=renamevars(orthgenes,'ortholog_ensg','id');
orthgenes=sortrows(orthgenes,'name');
orthgenes(orthgenes.id=="N/A",:)=[];