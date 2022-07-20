function [result, restab]= gene_tfidf(genes, counts, group, options)
arguments
    genes
    counts
    group
    options.N=50
    options.min_count=0.1
    options.min_in_freq=0.05
    options.max_out_freq=1
    options.do_hyge=1
    options.fdr_thr=0.05
    options.sortby="tfidf"
    options.sortdir="descend"
end

nCells=size(counts,2);

if isnumeric(group)
    group = "c"+string(group);
    group=categorical(group,natsort(unique(cellstr(group))));
end


if ~iscategorical(group)
    group=categorical(group);
end
group=removecats(group);
groupNames = categories(group);
nGroups=length(groupNames);

for i=1:nGroups
    if ~isnan(str2double(groupNames{i}))
        groupNames{i} = "c"+ groupNames{i};
    end
end

clusterSizes=countcats(group);

%grouped counts
[g, gn]=findgroups(group);
counts_thr=counts>=options.min_count;
groupCounts=splitapply(@(x)sum(x,2),counts_thr,g);
% groupCounts=splitapply(@(x)sum(x,2),counts,g);

keep=any(groupCounts,2);
groupCounts(~keep,:)=[];
genes(~keep,:)=[];

totalCounts=sum(groupCounts,2);

tf = groupCounts./clusterSizes;
ntf = (totalCounts-groupCounts)./(nCells-clusterSizes);
idf = log(nCells./totalCounts);
% idf = log((nCells-totalCounts)./totalCounts);

tfidf = tf.*idf;

%hypergeometric test for observed counts
if options.do_hyge
    p=ones(nnz(keep),nGroups);
    for i=1:length(groupNames)
        p(:,i)=hygecdf(groupCounts(:,i),nCells,totalCounts,clusterSizes(i),'upper');
    end
    p(p==0)=eps(0);
    
    [~,~,~,adj_p]=fdr_bh(p);
    sig=adj_p<options.fdr_thr;
    
    p(p>1)=1;
end

%sort genes by tfidf and retain top hits
%store the results in a table
restab=table();
for i=1:nGroups
        
    res=table;
    res.name=genes.name;
    res.in_freq=tf(:,i);
    res.out_freq=ntf(:,i);
    res.tot_freq=totalCounts./nCells;
    res.idf=idf;
    res.tfidf=tfidf(:,i);
    
    if options.do_hyge
        res.adj_p=adj_p(:,i);
        res(~sig(:,i),:)=[];
    end

    res(res.in_freq<options.min_in_freq,:)=[];
    res(res.out_freq>options.max_out_freq,:)=[];

    res=sortrows(res,'tfidf','descend');
    
    res=res(1:min(options.N,height(res)),:);

    res=sortrows(res,options.sortby,options.sortdir);

    result.(groupNames{i})=res;

    res.clust_id=repmat(string(groupNames{i}),height(res),1);
    restab=[restab;res];
end
% result=[];