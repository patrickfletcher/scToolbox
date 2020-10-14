function [dominant,specific,allgenes,dominantCondition, specificCondition]...
    =identifyDominantGenes(genes, selfin, otherin, P, params)

%self/other names must be same categories used to generate pairwise P values..

%TODO: check for & exclude internal dominance when self is a group? 
% - i.e. homogeneity constraint: strict co-dominance

%TODO: alternate interface - tcounts+categorical; compute the mean expr &
%prct expressing here. Also, compute P if not passed in...

%parameters
[self,other]=parse_input_pars(selfin,otherin);

nGenes=height(genes);

% whole group quantities, tests: P.anova, P.chi2
test_pAnova = P.anova < self.pthr;
test_pChi2 = P.chi2 < self.pprothr;

% group-wise quantities, tests: prct, expr
exprselfnames=strcat("expr_",self.names(:));
exprothernames=strcat("expr_",other.names(:));
prctselfnames=strcat("prct_",self.names(:));
prctothernames=strcat("prct_",other.names(:));

expr_self=genes{:,"expr_"+self.names(:)};
expr_other=genes{:,"expr_"+other.names(:)};
prct_self=genes{:,"prct_"+self.names(:)};
prct_other=genes{:,"prct_"+other.names(:)};

test_selfMinPrct=prct_self>self.minprct;
test_otherMaxPrct=all(prct_other<other.maxprct,2);

% group reductions
minSelfExpr=min(expr_self,[],2);
maxSelfExpr=max(expr_self,[],2);
minSelfPrct=min(prct_self,[],2);
maxSelfPrct=max(prct_self,[],2);
minOtherPrct=min(prct_other,[],2);
maxOtherPrct=max(prct_other,[],2);
minOtherExpr=min(expr_other,[],2);
maxOtherExpr=max(expr_other,[],2);

%make the following parameters:
exprESname='fc_expr';
prctESname='d_prct';

% pairwise quantities, tests: effect sizes, P values
fcCombs=allcomb(self.names,other.names); %nFcCombs=size(fcCombs,1);
pairwise_fc_expr=[];
pairwise_log2fc_expr=[];
pairwise_d_prct=[];
% pairwise_d_prct_rel=[];
pairwise_exprES=[]; 
pairwise_prctES=[]; 
pairwise_Pmc=[]; %store pairwise quantities in matrices for export
pairwise_Pz=[];
self_dominant=zeros(nGenes,length(self.names));
self_test_prctES=zeros(nGenes,length(self.names));
self_test_exprES=zeros(nGenes,length(self.names));

for i=1:length(self.names)
    
    %metrics
    this_fc_expr=expr_self(:,i)./expr_other; %one-centered, positive
    this_log2fc_expr=log2(this_fc_expr); %zero-centered - use log2 for easy double/half threhold
    this_d_prct=prct_self(:,i)-prct_other; %zero-centered
%     this_d_prct_rel=this_d_prct./maxOtherPrct;%zero-centered
    
    pairwise_fc_expr=[pairwise_fc_expr,this_fc_expr];
    pairwise_log2fc_expr=[pairwise_log2fc_expr,this_log2fc_expr];
    pairwise_d_prct=[pairwise_d_prct,this_d_prct];
%     pairwise_d_prct_rel=[pairwise_d_prct_rel,this_d_prct_rel];

    %effect sizes (TODO: add options to select different metrics)
%     switch params.prctMetric
%         case 'diff'
%         case 'reldiff'
%         case 'fc'
%         case 'log2fc'
%         case 'cohen_h'
%             h=2*asin(sqrt(prct_self/100)) - 2*asin(sqrt(prct_other/100));
%     end
%     switch params.exprMetric
%         case 'diff'
%         case 'reldiff'
%         case 'fc'
%         case 'log2fc'
%     end

    this_exprES=this_fc_expr;
    this_prctES=this_d_prct;
    
    pairwise_exprES=[pairwise_exprES,this_exprES];
    pairwise_prctES=[pairwise_prctES,this_prctES];
    
    %Pairwise p-values. loop only to get otherNameInCombo
    selfNameInCombo=any(strcmp(P.combs,self.names{i}),2);
    otherNameInCombo=any(ismember(P.combs,other.names),2);
    this_pmc=P.mc(:,selfNameInCombo&otherNameInCombo);
    this_pz=P.z(:,selfNameInCombo&otherNameInCombo);

    pairwise_Pmc=[pairwise_Pmc, this_pmc];
    pairwise_Pz=[pairwise_Pz, this_pz];

    test_thisExpr_P=test_pAnova & this_pmc<self.pthrpw;
    test_thisPrct_P=test_pChi2 & this_pz<self.pprothrpw;
    
%     switch params.Expr_ES_direction
%         case 'up'
            test_thisExpr_ES=this_exprES>self.fcExprThr;
            test_thisPrct_ES=this_prctES>self.prctESThr;
            test_pairwise_ESexpr_limit=this_exprES>self.minFcExpr; %lower bound for ES
            test_pairwise_ESprct_limit=this_prctES>self.minPrctEffect;
%         case 'down'
%             test_thisExpr_ES=this_exprES<self.fcExprThr;
%             test_thisPrct_ES=this_prctES<self.prctESThr;
%             test_pairwise_ESexpr_limit=this_exprES<self.maxFcExpr;
%             test_pairwise_ESprct_limit=this_prctES<self.maxPrctEffect;
%         case 'both'
%             test_thisExpr_ES=this_exprES>self.fcExprThr | this_exprES<self.fcExprThr;
%             test_thisPrct_ES=this_prctES>self.prctESThr | this_prctES<self.prctESThr;
%             test_pairwise_ESexpr_limit=this_exprES>self.minFcExpr | this_exprES<self.maxFcExpr;
%             test_pairwise_ESprct_limit=this_prctES>self.minPrctEffect | this_prctES<self.maxPrctEffect;
%     end
    
    
    % combine pairwise P value tests and effect-size tests: significant pairwise ES tests 
    test_thisExpr= test_thisExpr_P & test_thisExpr_ES;
    test_thisPrct= test_thisPrct_P & test_thisPrct_ES;
    
        
    %combine everything: apply minself, pairwise ES limits, and actual ES tests
    %for each test separately (report in case of interest):
    pairwise_Expr=test_selfMinPrct(:,i) & test_pairwise_ESprct_limit & test_pairwise_ESexpr_limit & test_thisExpr;
    pairwise_Prct=test_selfMinPrct(:,i) & test_pairwise_ESprct_limit & test_pairwise_ESexpr_limit & test_thisPrct;
    
    % how do we combine the two ES tests between this_self and each other?
%     switch params.pool_singlepairEStest_method
    switch self.combine_prct_expr
        case 'or'
            pairwise_combinedES=pairwise_Expr|pairwise_Prct; %permits either test
        case 'and'
            pairwise_combinedES=pairwise_Expr&pairwise_Prct; %strict
    end
    
    %NOW: combine all the pairwise tests. ALL by default. option for any?
    self_dominant(:,i)=all(pairwise_combinedES,2);
    self_test_exprES(:,i)=all(pairwise_Expr,2);
    self_test_prctES(:,i)=all(pairwise_Prct,2);
end

% gather extra info to report: to help understand what happened
% pariwise p-value reductions
max_Pmc=max(pairwise_Pmc,[],2);
max_Pz=max(pairwise_Pz,[],2);

% pariwise effect size reductions
min_fc_expr=min(pairwise_fc_expr,[],2);
max_fc_expr=max(pairwise_fc_expr,[],2);
min_log2fc_expr=min(pairwise_log2fc_expr,[],2);
max_log2fc_expr=max(pairwise_log2fc_expr,[],2);
min_d_prct=min(pairwise_d_prct,[],2);
max_d_prct=max(pairwise_d_prct,[],2);
max_abs_d_prct=max(abs(pairwise_d_prct),[],2);
% min_d_prct_rel=min(pairwise_d_prct_rel,[],2);
% max_d_prct_rel=max(pairwise_d_prct_rel,[],2);

%pairwise self pooling
% switch params.pool_selftest_method
if self.poolmethod=="all"
    pairwiseTestPooled=all(self_dominant,2);
    exprTestPooled=all(self_test_exprES,2);
    prctTestPooled=all(self_test_prctES,2);
else
    pairwiseTestPooled=any(self_dominant,2);
    exprTestPooled=any(self_test_exprES,2);
    prctTestPooled=any(self_test_prctES,2);
end

%assemble the final selections
dominantCondition = pairwiseTestPooled;
specificCondition = dominantCondition & test_otherMaxPrct;

%tables for output
allgenes=table();
allgenes.id=genes.id;
allgenes.name=genes.name;

allgenes.all_expr_test=exprTestPooled;
allgenes.all_prct_test=prctTestPooled;
allgenes.specific=specificCondition;

%effect sizes
allgenes.min_fc_expr=min_fc_expr;
allgenes.max_fc_expr=max_fc_expr;
allgenes.min_log2fc_expr=min_log2fc_expr;
allgenes.max_log2fc_expr=max_log2fc_expr;
allgenes.min_d_prct=min_d_prct;
allgenes.max_d_prct=max_d_prct;
% allgenes.min_d_prct_rel=min_d_prct_rel;
% allgenes.max_d_prct_rel=max_d_prct_rel;

%fc test
allgenes.p_anova=P.anova;
allgenes.max_pwp=max_Pmc;

%proportion test
allgenes.p_chi2=P.chi2;
allgenes.max_pwz=max_Pz;

%always output same format: if self size is 1, repeats...
% if length(self.names)>1
allgenes.min_self_expr=minSelfExpr;
allgenes.max_self_expr=maxSelfExpr;
allgenes.min_self_prct=minSelfPrct;
allgenes.max_self_prct=maxSelfPrct;
% else
%     allgenes.self_expr=minSelfExpr;
%     allgenes.self_prct=minSelfPrct;
% end

%specific (max other test)
allgenes.min_other_prct=minOtherPrct;
allgenes.max_other_prct=maxOtherPrct;
allgenes.min_other_expr=minOtherExpr;
allgenes.max_other_expr=maxOtherExpr;

allgenes.dominant=dominantCondition; %mainly for allgenes

%raw expression/prct values per type
allgenes=[allgenes,genes(:,...
    ismember(genes.Properties.VariableNames,[prctselfnames;prctothernames])...
   |ismember(genes.Properties.VariableNames,[exprselfnames;exprothernames]))];

%effect sizes
allgenes=[allgenes,array2table(pairwise_prctES,'variablenames',strcat(prctESname,'_',fcCombs(:,1),'_',fcCombs(:,2)))];
allgenes=[allgenes,array2table(pairwise_exprES,'variablenames',strcat(exprESname,'_',fcCombs(:,1),'_',fcCombs(:,2)))];

%pvals
allgenes=[allgenes,array2table(pairwise_Pmc,'variablenames',strcat('pmc_',fcCombs(:,1),'_',fcCombs(:,2)))];
allgenes=[allgenes,array2table(pairwise_Pz,'variablenames',strcat('pz_',fcCombs(:,1),'_',fcCombs(:,2)))];

dominant=allgenes(dominantCondition,:);
specific=allgenes(specificCondition,:);

% if length(self.names)>1
allgenes=allgenes(~(allgenes.max_self_expr==0)&~(allgenes.max_other_expr==0),:);
% else
% allgenes=allgenes(~(allgenes.self_expr==0)&~(allgenes.max_other_expr==0),:);
% end

end


function [selfpars,otherpars]=parse_input_pars(selfin,otherin)
selfpars.names={'self'};
selfpars.minprct=0;
selfpars.prctMetric='diff';
selfpars.prctESThr=30;
selfpars.minPrctEffect=0;
selfpars.fcExprThr=0;
selfpars.minFcExpr=0;
selfpars.combine_prct_expr='or';
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

end

