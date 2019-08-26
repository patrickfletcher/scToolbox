function geneTable=getExpression(genes,ncounts,tcounts,factor1,factor2)
%get gene expression for groups defined by two levels of factors

factor1Names=categories(factor1);

do2factor=false;
if exist('factor2','var')&&~isempty(factor2)
    factor2Names=categories(factor2);
    do2factor=true;
end

geneTable=genes(:,{'id','name'});
for i=1:length(factor1Names)
    if do2factor
        for j=1:length(factor2Names)
            thisName=[factor1Names{i},'_',factor2Names{j}];
            thisGroup=factor1==factor1Names(i) & factor2==factor2Names{j};
            geneTable.(['prct_',thisName])=sum(tcounts(:,thisGroup)>genes.thr,2)./nnz(thisGroup)*100;
            geneTable.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
        end
    else
        thisName=factor1Names{i};
        thisGroup=factor1==factor1Names(i);
        geneTable.(['prct_',thisName])=sum(tcounts(:,thisGroup)>genes.thr,2)./nnz(thisGroup)*100;
        geneTable.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
    end
end

%reorder the columns in alphabetical order
geneTable=geneTable(:,sort(geneTable.Properties.VariableNames));