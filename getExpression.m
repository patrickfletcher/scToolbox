function geneTable=getExpression(genes,ncounts,tcounts,factor1,factor2,threshgroup)
%get gene expression for groups defined by two levels of factors
%
% assumes tcounts is already threshold-subtracted, unless a threshold-group
% is passed in...
factor1Names=categories(factor1);

do2factor=false;
if exist('factor2','var')&&~isempty(factor2)
    factor2Names=categories(factor2);
    do2factor=true;
end

doThreshold=false;
if exist('threshgroup','var')&&~isempty(threshgroup)
    threshgroup=categorical(threshgroup);
    threshgroups=categories(threshgroup);
    doThreshold=true;
end
if doThreshold
    if length(threshgroups)>1
        for i=1:length(threshgroups)
            thisgroup=threshgroup==threshgroups{i};
            tcounts(:,thisgroup)=tcounts(:,thisgroup)-genes.("thr_"+threshgroups{i});
        end
    else
        tcounts=tcounts-genes.thr;
    end
end

geneTable=genes(:,{'id','name'});
for i=1:length(factor1Names)
    if do2factor
        for j=1:length(factor2Names)
            thisName=[factor1Names{i},'_',factor2Names{j}];
            thisGroup=factor1==factor1Names(i) & factor2==factor2Names{j};
            geneTable.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
            geneTable.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
        end
    else
        thisName=factor1Names{i};
        thisGroup=factor1==factor1Names(i);
        geneTable.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
        geneTable.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
    end
end

%reorder the columns in alphabetical order
geneTable=geneTable(:,sort(geneTable.Properties.VariableNames));