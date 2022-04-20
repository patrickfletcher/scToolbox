function [geneTable, E, P, factorNames]=getExpression(genes,ncounts,tcounts,factor,options)
arguments
    genes
    ncounts
    tcounts
    factor
    options.factor2 = []
    options.threshgroup = []
    options.method = 'mean'
    options.trim_prct = 0
    options.doPrct = true
    options.doPooled = false %true: Cat1, nonCat1, Cat2, ... 
    options.only_expressing = false;
end
%[geneTable, E, P, factorNames]=getExpression(genes,ncounts,tcounts,factor,options)
%get gene expression for groups defined by two levels of factors
%
% assumes tcounts is already threshold-subtracted, unless a threshold-group
% is passed in...

%TODO: optional ncounts/tcounts input for prct (i.e. part of thresholding)

%remove cats first?
if ~iscategorical(factor)
    factor=categorical(factor);
end
factor=removecats(factor);

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

switch options.method
    case 'trimean'
        expr_fun=@(x) mean(prctile(x,[25,50,50,75],2),2); %fastest
    case 'trimmean'
        expr_fun=@(x) trimmean(x, options.trim_prct, 2);
    case 'sum'
        expr_fun=@(x) sum(x,2);
    otherwise
        expr_fun=eval("@(x)"+options.method+"(x,2,'omitnan')");
end

numCells=size(ncounts,2);
E=zeros(size(tcounts,1),length(factorNames));
P=zeros(size(tcounts,1),length(factorNames));
if options.doPooled
    Eother=zeros(size(tcounts,1),length(factorNames));
    Pother=zeros(size(tcounts,1),length(factorNames));
end
for i=1:length(factorNames)
    thisGroup=factor==factorNames(i);
    
    N1=ncounts(:,thisGroup);
    if options.only_expressing
        N1(N1==0)=nan;
    end

    E(:,i) = expr_fun(N1);
    
    if options.doPooled
        N2=ncounts(:,~thisGroup);
        Eother(:,i)=expr_fun(N2);
    end
    
    if options.doPrct
        T1=tcounts(:,thisGroup);
        T1=T1>0;
        c1=nnz(thisGroup);
        P(:,i)=sum(T1,2)./c1*100;
%         if options.doPooled
%             T2=tcounts(:,~thisGroup);
%             T2=T2>0;
%             c2=numCells-c1;
%             Pother(:,i)=sum(T2,2)./c2*100;
%         end
    end
end

geneTable=genes(:,{'id','name'});
geneTable=[geneTable,array2table(E,'VariableNames',"expr_"+string(factorNames))];
if options.doPooled
    geneTable=[geneTable,array2table(Eother,'VariableNames',"expr_not_"+string(factorNames))];
end
if options.doPrct
    geneTable=[geneTable,array2table(P,'VariableNames',"prct_"+string(factorNames))];
%     if options.doPooled
%         geneTable=[geneTable,array2table(Pother,'VariableNames',"prct_not_"+string(factorNames))];
%     end
end

end

% function TM = trimean(x, dim)
% arguments
%     x
%     dim = 1
% end
% med=median(x,dim);
% qts=prctile(x,[25,75],dim);
% TM=sum(2*med+qts, dim)/4;
% end

    
%     switch options.method
%         case 'mean'
%             thisExpr=mean(ncounts(:,thisGroup),2);
%         case 'median'
%             thisExpr=median(ncounts(:,thisGroup),2);
%         case 'geomean'
%             thisExpr=geomean(ncounts(:,thisGroup),2);
%         case 'harmmean'
%             thisExpr=harmmean(ncounts(:,thisGroup),2);
%         case 'trimean'
%             thisExpr=mean(prctile(ncounts(:,thisGroup),[25,50,50,75],2),2);
%         case 'trimmean'
%             thisExpr=trimmean(ncounts(:,thisGroup),options.trim_prct,2);
%     end
%     E(:,i) = thisExpr;

% %with grouped summary - very slow for large matrices...
% %average expression
% switch options.method
%     case 'trimean'
%         expr_fun=@(x) mean(prctile(x,[25,50,50,75]));
%     case 'trimmean'
%         expr_fun=@(x) trimmean(x, options.trim_prct);
%     otherwise
%         expr_fun=options.method;
% end
% [E,factorNames]=groupsummary(ncounts',factor,expr_fun);
% E=E';
% factorNames=factorNames';
% 
% % %percent expressing
% P=groupsummary(tcounts'>0,factor,'sum');
% P=P.*100./groupcounts(factor);
% P=P';


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
