function geneTable=getExpression(genes,ncounts,tcounts,factor1,factor2,threshgroup)
%geneTable=getExpression(genes,ncounts,tcounts,factor1,factor2,threshgroup)
%
%get gene expression for groups defined by two levels of factors
%
% assumes tcounts is already threshold-subtracted, unless a threshold-group
% is passed in...

%TODO: use stratify_factors

%remove cats first?
if ~iscategorical(factor1)
    factor1=categorical(factor1);
end
factor1=removecats(factor1);
factor1Names=categories(factor1);

do2factor=false;
if exist('factor2','var')&&~isempty(factor2)
    %remove cats first?
    if ~iscategorical(factor2)
        factor2=categorical(factor2);
    end
    factor2=removecats(factor2);
    factor2Names=categories(factor2);
    do2factor=true;
end

doThreshold=false;
if exist('threshgroup','var')
    if isscalar(threshgroup)
        if threshgroup==1 %will use genes.thr
            threshgroups=1;
            doThreshold=true;
        end
    elseif length(threshgroup)==size(tcounts,2)
        threshgroup=categorical(threshgroup);
        threshgroups=categories(threshgroup);
        doThreshold=true;
    end
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

geneExpr=table();
genePrct=table();
for i=1:length(factor1Names)
    if do2factor
        for j=1:length(factor2Names)
            thisName=[factor1Names{i},'_',factor2Names{j}];
            thisGroup=factor1==factor1Names(i) & factor2==factor2Names{j};
            genePrct.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
            geneExpr.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
        end
    else
        thisName=factor1Names{i};
        thisGroup=factor1==factor1Names(i);
        genePrct.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
        geneExpr.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
    end
end

geneTable=[genes(:,{'id','name'}),geneExpr,genePrct];