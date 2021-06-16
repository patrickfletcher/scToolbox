function geneTable=getExpression(genes,ncounts,tcounts,factor,options)
arguments
    genes
    ncounts
    tcounts
    factor
    options.factor2 = []
    options.threshgroup = []
    options.method = 'mean'
    options.trim_prct = 0
end
%get gene expression for groups defined by two levels of factors
%
% assumes tcounts is already threshold-subtracted, unless a threshold-group
% is passed in...

%TODO: remove second factor?
%TODO: option for either mean/prct?
%TODO: use accumarray/groupcounts

%remove cats first?
if ~iscategorical(factor)
    factor=categorical(factor);
end
factor=removecats(factor);
factor1Names=categories(factor);

do2factor=false;
if ~isempty(options.factor2)
    %remove cats first?
    factor2=categorical(options.factor2);
    factor2=removecats(factor2);
    do2factor=true;
end

doThreshold=false;
if ~isempty(options.threshgroup)
    if isscalar(options.threshgroup)
        if options.threshgroup==1 %will use genes.thr
            threshgroups=1;
            doThreshold=true;
        end
    elseif length(options.threshgroup)==size(tcounts,2) %should really compute thr here, not expect a column in genes
        threshgroup=categorical(options.threshgroup);
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


if do2factor
    factor = stratify_factors(factor, factor2);
    factor = removecats(factor);
end
factorNames = categories(factor);

geneExpr=table();
genePrct=table();
for i=1:length(factorNames)
    thisName=factorNames{i};
    thisGroup=factor==factorNames(i);
    genePrct.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
    switch options.method
        case 'mean'
            thisExpr=mean(ncounts(:,thisGroup),2);
        case 'median'
            thisExpr=median(ncounts(:,thisGroup),2);
        case 'geomean'
            thisExpr=geomean(ncounts(:,thisGroup),2);
        case 'harmmean'
            thisExpr=harmmean(ncounts(:,thisGroup),2);
        case 'trimean'
            thisExpr=mean(prctile(ncounts(:,thisGroup),[25,50,50,75],2),2);
        case 'trimmean'
            thisExpr=trimmean(ncounts(:,thisGroup),options.trim_prct,2);
    end
    geneExpr.(['expr_',thisName]) = thisExpr;
end

geneTable=[genes(:,{'id','name'}),geneExpr,genePrct];


% for i=1:length(factor1Names)
%     if do2factor
%         for j=1:length(factor2Names)
%             thisName=[factor1Names{i},'_',factor2Names{j}];
%             thisGroup=factor1==factor1Names(i) & factor2==factor2Names{j};
%             genePrct.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
%             geneExpr.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
%         end
%     else
%         thisName=factor1Names{i};
%         thisGroup=factor1==factor1Names(i);
%         genePrct.(['prct_',thisName])=sum(tcounts(:,thisGroup)>0,2)./nnz(thisGroup)*100;
%         switch method
%             case 'mean'
%         geneExpr.(['expr_',thisName])=mean(ncounts(:,thisGroup),2);
%             case 'median'
%             case 'geomean'
%             case 'harmmean'
%             case 'trimean'
%             case 'trimmean'
%         end
%     end
% end
