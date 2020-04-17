function [dominant,specific, dominantCondition, specificCondition]...
    =identifyDominantGenes(genes, selfin, otherin, P)
%
%self/other names must be categories used to generate pairwise P values..

%put combo names as varnames in P table? eg. E_Le -> get tokens surrounding
%"_" and then search for each self.name.

%if doing self vs. other (one comparison):
% - P should be a column vector
% - "genes" should include prct_self/prct_other/expr_...
% - combinations not needed - set it to {'self','other'}?

%BUG: sometimes fctest 

%TODO: internally convert foldchange to log10foldchange?

%TODO: check for & exclude internal dominance when self is a group? i.e. homogeneity constraint

doProportionTest=false;
if isfield(P,'chi2')
    doProportionTest=true;
end

%parameters
selfpars.names={'self'};
selfpars.minprct=0;
selfpars.prctMetric='diff';
selfpars.prctESThr=30;
selfpars.minPrctEffect=0;
selfpars.fcExprThr=0;
selfpars.minFcExpr=0;
selfpars.prctOrExpr=false;
selfpars.pthr=0.1;
selfpars.pthrpw=0.1;
selfpars.pprothr=0.1;
selfpars.pprothrpw=0.1;
selfpars.poolmethod={'all'};

otherpars.names={'other'};
otherpars.maxprct=Inf;
otherpars.poolmethod={'all'};

%overwrite field values with any that were passed in
for f=intersect(fieldnames(selfpars),fieldnames(selfin))'
    selfpars.(f{1})=selfin.(f{1});
end
for f=intersect(fieldnames(otherpars),fieldnames(otherin))'
    otherpars.(f{1})=otherin.(f{1});
end

% if any(~ismember(fieldnames(selfin),fieldnames(self)))
%     warning('unused fields in self-struct')
% end
% if any(~ismember(fieldnames(otherin),fieldnames(other)))
%     warning('unused fields in self-struct')
% end

nGenes=height(genes);

%extract the relevant prct, expr, and P values
prctselfnames=strcat('prct_',selfpars.names(:));
exprselfnames=strcat('expr_',selfpars.names(:));
prctothernames=strcat('prct_',otherpars.names(:));
exprothernames=strcat('expr_',otherpars.names(:));
selfPrct=genes{:,prctselfnames};
otherPrct=genes{:,prctothernames};
selfExpr=genes{:,exprselfnames};
otherExpr=genes{:,exprothernames};

%all possible combinations of self with other
fcCombs=allcomb(selfpars.names,otherpars.names);
nFcCombs=size(fcCombs,1);

%get p-values for comparisons of interest
thisPmc=ones(nGenes,nFcCombs);
thisPz=zeros(nGenes,nFcCombs); %defaults to below thresh if ~doProportionTest
ix=1;
for i=1:length(selfpars.names)
    selfNameInCombo=any(strcmp(P.combs,selfpars.names{i}),2);
    for j=1:length(otherpars.names)
        otherNameInCombo=any(strcmp(P.combs,otherpars.names{j}),2);
        thisPmc(:,ix)=P.mc(:,selfNameInCombo&otherNameInCombo);
        if doProportionTest
            thisPz(:,ix)=P.z(:,selfNameInCombo&otherNameInCombo);
        end
        ix=ix+1; %(i-1)*length(otherpars.names)+j
    end
end

%apply expr and prct test pairwise.  All pairs containing a single "self" group must pass. pool self groups using "all"
%or "any" as needed
dPrct=zeros(nGenes,nFcCombs);
dExpr=ones(nGenes,nFcCombs);
fcPrct=ones(nGenes,nFcCombs);
fcExpr=ones(nGenes,nFcCombs);
% pairwiseTests=false(nGenes,nFcCombs);
selfDominant=false(nGenes,length(selfpars.names));
selfExprTest=false(nGenes,length(selfpars.names));
selfPrctTest=false(nGenes,length(selfpars.names));
% ix=1;
for i=1:length(selfpars.names)
    pairwiseTests=false(nGenes,length(otherpars.names));
    pairwiseExprTests=false(nGenes,length(otherpars.names));
    pairwisePrctTests=false(nGenes,length(otherpars.names));
    for j=1:length(otherpars.names)
        
        ix=(i-1)*length(otherpars.names)+j;
        pmc=thisPmc(:,ix);
        pz=thisPz(:,ix);
        
        prct1=selfPrct(:,i);
        prct2=otherPrct(:,j);
        expr1=selfExpr(:,i);
        expr2=otherExpr(:,j);
        
        dExpr(:,ix)=expr1-expr2;
        
        [pairwiseTests(:,j),dPrct(:,ix),fcPrct(:,ix),fcExpr(:,ix),pairwiseExprTests(:,j),pairwisePrctTests(:,j)]...
            =compareTwoGroups(P.anova,pmc,P.chi2,pz,prct1,prct2,expr1,expr2,selfpars);
        
%         ix=ix+1;
    end
    
    %thresholds on self%, minimum effect sizes
    selfMinTest(:,i)=prct1>selfpars.minprct;
    
    thisDprct=selfPrct(:,i)-otherPrct;
    pairwisePrctEffectTests=thisDprct>selfpars.minPrctEffect;
    
    thisFCexpr=selfExpr(:,i)./otherExpr;
    pairwiseFCexprTests=thisFCexpr>selfpars.minFcExpr;
    
    pairwiseDominance=selfMinTest(:,i) & pairwisePrctEffectTests & pairwiseFCexprTests & (pairwiseExprTests|pairwisePrctTests);
    
    %pairwise other pooling - any reason to use "any"?
    selfDominant(:,i)=all(pairwiseDominance,2);
%     selfDominant(:,i)=all(pairwiseTests,2);
    minPrctEffectTest(:,i)=all(pairwisePrctEffectTests,2);
    minFCexprTest(:,i)=all(pairwiseFCexprTests,2);
    selfExprTest(:,i)=all(pairwiseExprTests,2);
    selfPrctTest(:,i)=all(pairwisePrctTests,2);
    
end

maxPairwiseP=max(thisPmc,[],2);
minDprct=min(dPrct,[],2);
minDexpr=min(dExpr,[],2);
minFCprct=min(fcPrct,[],2);
minFCexpr=min(fcExpr,[],2);
minSelfPrct=min(selfPrct,[],2);
minSelfExpr=min(selfExpr,[],2);
maxOtherPrct=max(otherPrct,[],2);

%pairwise self pooling
if selfpars.poolmethod=="all"
    pairwiseTestPooled=all(selfDominant,2);
    exprTestPooled=all(selfExprTest,2);
    prctTestPooled=all(selfPrctTest,2);
else
    pairwiseTestPooled=any(selfDominant,2);
    exprTestPooled=any(selfExprTest,2);
    prctTestPooled=any(selfPrctTest,2);
end

% %ANOVA-type tests
% pCriterionANOVA = P.anova<selfpars.pthr;
if doProportionTest
    %     pCriterionANOVA=pCriterionANOVA & P.chi2<selfpars.pthr;
    maxPairwisePz=max(thisPz,[],2);
end

%other percent tests
otherMaxPrctCriterionPool=maxOtherPrct<=otherpars.maxprct;

%assemble the final selections
dominantCondition = pairwiseTestPooled;
% dominantCondition = pCriterionANOVA & pairwiseTestPooled;
specificCondition = dominantCondition & otherMaxPrctCriterionPool;

%tables for output
dominant=table();
dominant.id=genes.id(dominantCondition);
dominant.name=genes.name(dominantCondition);
if length(selfpars.names)>1
    dominant.min_self_prct=minSelfPrct(dominantCondition);
    dominant.min_self_expr=minSelfExpr(dominantCondition);
else
    dominant.self_prct=minSelfPrct(dominantCondition);
    dominant.self_expr=minSelfExpr(dominantCondition);
end
%specific (max other test)
dominant.max_other_prct=maxOtherPrct(dominantCondition);

%effect sizes
dominant.min_fc_expr=minFCexpr(dominantCondition);
dominant.min_d_prct=minDprct(dominantCondition);

%fc test
dominant.all_fc_test=exprTestPooled(dominantCondition);
dominant.p_anova=P.anova(dominantCondition);
dominant.max_pwp=maxPairwiseP(dominantCondition);

%proportion test
dominant.all_prct_test=prctTestPooled(dominantCondition);
if doProportionTest
    dominant.p_chi2=P.chi2(dominantCondition);
    dominant.max_pwz=maxPairwisePz(dominantCondition);
end

%NOTE: all_fc_test=0 & all_prct_test=0 means for some pairs, one of the two
%tests passed but not the other, so that not ALL pairs passed FC and not
%ALL pairs passed PRCT.

% dominant.min_d_expr=minDexpr(dominantCondition);
% dominant.min_fc_prct=minFCprct(dominantCondition);

%raw expression/prct values per type
dominant=[dominant,genes(dominantCondition,...
    ismember(genes.Properties.VariableNames,[prctselfnames;prctothernames])...
    |ismember(genes.Properties.VariableNames,[exprselfnames;exprothernames]))];

%effect sizes
dominant=[dominant,array2table(dPrct(dominantCondition,:),'variablenames',strcat('dprct_',fcCombs(:,1),'_',fcCombs(:,2)))];
% dominant=[dominant,array2table(fcPrct(dominantCondition,:),'variablenames',strcat('fcprct_',fcCombs(:,1),'_',fcCombs(:,2)))];
dominant=[dominant,array2table(fcExpr(dominantCondition,:),'variablenames',strcat('fcexpr_',fcCombs(:,1),'_',fcCombs(:,2)))];

%pvals
dominant=[dominant,array2table(thisPmc(dominantCondition,:),'variablenames',strcat('pmc_',fcCombs(:,1),'_',fcCombs(:,2)))];
if doProportionTest
    dominant=[dominant,array2table(thisPz(dominantCondition,:),'variablenames',strcat('pz_',fcCombs(:,1),'_',fcCombs(:,2)))];
end

specific=dominant(specificCondition(dominantCondition),:);

end


function [pwtest,dprct,fcprct,fcexpr,exprTest,prctTest]=compareTwoGroups(Panova,pmc,Pchi2,pz,prct1,prct2,expr1,expr2,par)
dprct=prct1-prct2;
fcprct=prct1./prct2;  %fcprct(prct1==0 & prct2==0)=1; %ie. no fold change: 0==0  - nan's don't compare anyways..
fcexpr=expr1./expr2;  %fcexpr(expr1==0 & expr1==0)=1;

switch par.prctMetric
    case 'diff'
        prctEffectSize=dprct;
    case 'fc'
        prctEffectSize=fcprct;
end

% exprTest=fcexpr>=par.fcExprThr;
exprTest=Panova<par.pthr & pmc<par.pthrpw & fcexpr>=par.fcExprThr;

% prctTest=prctEffectSize>=par.prctESThr;
prctTest=Pchi2<par.pprothr & pz<par.pprothrpw & prctEffectSize>=par.prctESThr;

if par.prctOrExpr
    pwtest=exprTest | prctTest;
else
    pwtest=exprTest & prctTest;
end

pwtest = pwtest & prct1>=par.minprct & fcexpr>par.minFcExpr & prctEffectSize>par.minPrctEffect;

%     if par.prctOrExpr
%         sufficientEffectSize=prctEffectSize>=par.prctESThr | fcexpr>=par.fcExprThr;
%     else
%         sufficientEffectSize=prctEffectSize>=par.prctESThr & fcexpr>=par.fcExprThr;
%     end
%
%     pwtest = prct1>=par.minprct & pmc<=par.pthr & pz<=par.pthr & sufficientEffectSize;
end


% %pvalue tests
% pCriterionANOVA = Panova<self.pthr;
% pCriterion = thisP<self.pthrpw;
% pCriterionPool = all(pCriterion,2);
%
% if doPZ
% pZCriterion = thisPz<self.pthr;
% pZCriterionPool = all(pZCriterion,2);
% end
%
% %fold change tests
% % selfDPrctCriterion= dPrct >= self.dPrctThr;
% selfFCprctCriterion= fcPrct >= self.fcPrctThr;
% selfFCexprCriterion= fcExpr >= self.fcExprThr;
% selfFCprctCriterionPool=all( fcPrct >= self.fcPrctThr, 2);
% selfFCexprCriterionPool=all( fcExpr >= self.fcExprThr, 2);
% if self.prctOrExpr
%     selfFCCriterion=selfFCprctCriterion | selfFCexprCriterion;
% else
%     selfFCCriterion=selfFCprctCriterion & selfFCexprCriterion;
% end

% %self percent tests
% minSelfPrct=min(selfPrct,[],2);
% selfMinPrctCriterionPool=minSelfPrct>=self.minprct;
%
% selfMinPrctCriterionRep=false(size(selfFCCriterion));
% selfMinPrctCriterion=selfPrct>self.minprct;
% for i=1:nFcCombs
%     selfIdx=strcmp(self.names,fcCombs{i,1});
%     selfMinPrctCriterionRep(:,i)=selfMinPrctCriterion(:,selfIdx);  %replicate the self % test for each pairwise test
% end
% selfMinPrctCriterion=selfPrct>=self.minprct;
% selfMinPrctCriterion=all(selfPrct>=self.minprct,2)

%self homogeneity if self is a group
% if length(self.names)>1
%     selfCombs=combnk(self.names,2);
%     nCombs=size(selfCombs,1);
%     for i=1:nCombs
%         selfIdx=strcmp(self.names,fcCombs{i,1});
%         otherIdx=strcmp(other.names,fcCombs{i,2});
%         selfDiff(:,i)=selfPrct(:,selfIdx)-otherPrct(:,otherIdx);
%     end
% end


%collect self tests - this way, each pairwise test is also coupled with the self min% test, so that genes with <self%
%but passing FC test are not included.
% if doPZ
% pairwiseTests=pCriterion & pZCriterion & selfFCCriterion & selfMinPrctCriterionRep;
% else
% pairwiseTests=pCriterionANOVA &pCriterion & selfFCCriterion & selfMinPrctCriterionRep;
% end
