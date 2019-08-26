function geneTable=getExpressionBySex(genes,ncounts,tcounts,factor1,sex)
%get gene expression for groups defined by two levels of factors

factor1Names=categories(factor1);

geneTable=genes(:,{'id','name'});
for i=1:length(factor1Names)
    thisName=[factor1Names{i},'_F'];
    thisGroup=factor1==factor1Names(i) & sex=="F";
    geneTable.(['prct_',thisName])=sum(tcounts(:,thisGroup)>genes.thr_F,2)./nnz(thisGroup)*100;
    geneTable.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
        
    thisName=[factor1Names{i},'_M'];
    thisGroup=factor1==factor1Names(i) & sex=="M";
    geneTable.(['prct_',thisName])=sum(tcounts(:,thisGroup)>genes.thr_M,2)./nnz(thisGroup)*100;
    geneTable.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
end

%reorder the columns in alphabetical order
geneTable=geneTable(:,sort(geneTable.Properties.VariableNames));