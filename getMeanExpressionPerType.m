function geneData=getMeanExpressionPerType(X,genes,group,nonzeroflag)
%returns a table with vars: gene, thr, %expressing per group, log10(mean(ncounts)+1) per group

if ~exist('nonzeroflag','var')
    nonzeroflag=false;
end

geneData=table();
geneData.gene=genes(:);

if ~iscategorical(group)
    group=categorical(group);
end
groupNames=categories(group);

if nonzeroflag
    X(X==0)=nan;
end

meanExpr=zeros(length(genes),length(groupNames));
for i=1:length(groupNames)
    meanExpr(:,i)=mean(X(:,group==groupNames(i)),2,'omitnan');
end

%optional?
% meanExpr=log10(meanExpr+1);

geneData=[geneData,array2table(meanExpr,'VariableNames',strcat('expr_',groupNames))];