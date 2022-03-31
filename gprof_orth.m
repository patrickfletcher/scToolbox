function [orthgenes, orthTab] = gprof_orth(query, source, target)

args=pyargs('base_url','https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15/');
gp=py.gprofiler.GProfiler(args);

pyquery=py.list(cellstr(query(:)'));
gOrthArgs=pyargs('organism',source,'query',pyquery,'target',target);

res=gp.orth(gOrthArgs); %returns list of dicts
res=Core_py2matlab(res); 
res=[res{:}];
orthTab=struct2table(res);

[orthgenes,ia]=unique(orthTab(:,{'ortholog_ensg','name'}),'rows','stable');
% orthTab.ia=ia;
orthgenes=renamevars(orthgenes,'ortholog_ensg','id');
orthgenes=sortrows(orthgenes,'name');
orthgenes(orthgenes.id=="N/A",:)=[];

if iscell(orthTab.description)
dclass=cellfun(@(x)class(x),orthTab.description,'UniformOutput',false);
orthTab.description(ismember(dclass,'py.NoneType'))={''};
end
orthTab=standardizeMissing(orthTab,{'N/A','nan'});