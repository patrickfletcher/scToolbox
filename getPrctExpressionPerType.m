function geneData=getPrctExpressionPerType(tcounts,gene_list,group,thresh)
%threshold and tcounts must match

geneData=table();
geneData.gene=gene_list(:);

%thresh is zero by default (no arg passed). Could be a scalar, or list of genewise thresholds. Can also pass char
%thresh method name (eg otsu)
if ~exist('thresh','var')
    thresh=zeros(size(gene_list(:)));
elseif isnumeric(thresh) && isscalar(thresh)
    thresh=thresh*ones(size(gene_list(:)));
elseif ischar(thresh)
    switch lower(thresh)
        case 'otsu'
            thresh=geneThresholdOtsu(tcounts);
    end
end

geneData.thr=thresh(:);

if ~iscategorical(group)
    group=categorical(group);
end
group=removecats(group);
groupNames=categories(group);
groupCounts=countcats(group);

prctPerType=zeros(length(gene_list),length(groupNames));
for i=1:length(groupNames)
    prctPerType(:,i)=sum(tcounts(:,group==groupNames(i))>thresh,2)./groupCounts(i)*100;
end

geneData=[geneData,array2table(prctPerType,'VariableNames',strcat('prct_',groupNames))];