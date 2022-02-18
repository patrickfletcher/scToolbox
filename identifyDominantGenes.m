function [dominant,specific,allgenes,dominantCondition, specificCondition]...
    =identifyDominantGenes(genes, selfin, otherin, P, paramsin)

%self/other names must be same categories used to generate pairwise P values..

% if P not supplied, just check the effect sizes!

%TODO: check for & exclude internal dominance when self is a group? 
% - i.e. homogeneity constraint: strict co-dominance

%TODO: alternate interface - tcounts+categorical; compute the mean expr &
%prct expressing here. Also, compute P if not passed in...

%TODO: toggle which ES to use (or combo)

%TODO: use arguments block.  Multiple structs can't have same fieldnames
% selfpars.names={'self'};
% selfpars.minprct=0;
% selfpars.maxprct=Inf;
% selfpars.poolmethod='all';
% selfpars.some_N=1; %as in for "any"
% otherpars.names={'other'};
% otherpars.minprct=0;
% otherpars.maxprct=Inf;
% otherpars.poolmethod='all';
% otherpars.some_N=1; %as in for "any"
% params.direction='any';
% params.limitMode='pairwise';
% params.combine_prct_expr='or';
% params.specificMode='uniform'; %or pairwise
% params.prctMetric='diff';
% params.prctESThr=20;
% params.fcExprThr=1;
% params.minPrctEffect=-100;
% params.minFcEffect=0;
% params.pthr=0.1;
% params.pthrpw=0.1;
% params.pprothr=0.1;
% params.pprothrpw=0.1;
% params.sortval='none';
% params.sortdir='descend';

%parameters
[self, other, params]=parse_input_pars(selfin,otherin, paramsin);

nGenes=height(genes);

% whole group quantities, tests: P.anova, P.chi2
Panova=zeros(nGenes,1);
Pchi2=zeros(nGenes,1);

%TODO: switch is based on presence of anova/chi2. probably a better way
usePexpr=false;
if ~isempty(P) && isfield(P,'anova')
    Panova = P.anova;
    usePexpr=true;
end
usePprop=false;
if ~isempty(P) && isfield(P,'chi2')
    Pchi2 = P.chi2;
    usePprop=true;
end
    
test_pAnova=Panova < params.pthr;
test_pChi2=Pchi2 < params.pprothr;

% group-wise quantities, tests: prct, expr
exprselfnames=strcat("expr_",self.names(:));
exprothernames=strcat("expr_",other.names(:));
prctselfnames=strcat("prct_",self.names(:));
prctothernames=strcat("prct_",other.names(:));

expr_self=genes{:,"expr_"+self.names(:)};
expr_other=genes{:,"expr_"+other.names(:)};
prct_self=genes{:,"prct_"+self.names(:)};
prct_other=genes{:,"prct_"+other.names(:)};

test_selfMinPrct=prct_self>self.minprct; %for up
test_otherMaxPrct=prct_other<other.maxprct;  %for up SPEC: all vs pairwise?
test_selfMaxPrct=prct_self<self.maxprct; %for down
test_otherMinPrct=prct_other>other.minprct;

%bug with self.poolmethod='any'/'some': need to apply the minPrct & specific tests
%per self.name in that case
switch params.direction
    case 'up'
        test_minPrct=all(test_selfMinPrct,2);
        test_specific=all(test_otherMaxPrct,2);
    case 'down'
        test_minPrct=all(test_otherMinPrct,2);
        test_specific=all(test_selfMaxPrct,2);
    case 'both'
        test_minPrct=all(test_selfMinPrct,2)|all(test_otherMinPrct,2);
        test_specific=all(test_otherMaxPrct,2)|all(test_selfMaxPrct,2);
end

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
pairwise_Pmc=[]; 
pairwise_Pz=[];
dominant_otherpooled=zeros(nGenes,length(self.names));
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
    % really should explore different % effect size metrics
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
    
    %Pairwise p-values. in loop because need otherNameInCombo
    if usePexpr||usePprop
        selfNameInCombo=any(strcmp(P.combs,self.names{i}),2);
        otherNameInCombo=any(ismember(P.combs,other.names),2);
    end
    this_pmc=zeros(size(this_exprES)); %makes P-val tests always true
    if usePexpr
        this_pmc=P.mc(:,selfNameInCombo&otherNameInCombo);
    end
    this_pz=zeros(size(this_prctES));
    if usePprop
        this_pz=P.z(:,selfNameInCombo&otherNameInCombo);
    end
    
    pairwise_Pmc=[pairwise_Pmc, this_pmc];
    pairwise_Pz=[pairwise_Pz, this_pz];

    test_thisExpr_P=test_pAnova & this_pmc<params.pthrpw;
    test_thisPrct_P=test_pChi2 & this_pz<params.pprothrpw;
    
    switch params.direction
        case 'up'
            test_thisExpr_ES=this_exprES>params.fcExprThr;
            test_thisPrct_ES=this_prctES>params.prctESThr;
            test_pairwise_ESexpr_limit=this_exprES>params.minFcEffect; %lower bound for ES
            test_pairwise_ESprct_limit=this_prctES>params.minPrctEffect;
        case 'down'
            test_thisExpr_ES=this_exprES<1/params.fcExprThr;
            test_thisPrct_ES=this_prctES<-params.prctESThr;
            test_pairwise_ESexpr_limit=this_exprES<1/params.minFcEffect;
            test_pairwise_ESprct_limit=this_prctES<-params.minPrctEffect;
        case 'both'
            test_thisExpr_ES=this_exprES>params.fcExprThr | this_exprES<1/params.fcExprThr;
            test_thisPrct_ES=this_prctES>params.prctESThr | this_prctES<-params.prctESThr;
            test_pairwise_ESexpr_limit=this_exprES>params.minFcEffect | this_exprES<1/params.minFcEffect;
            test_pairwise_ESprct_limit=this_prctES>params.minPrctEffect | this_prctES<-params.minPrctEffect;
    end

    % combine pairwise P value tests and effect-size tests: significant pairwise ES tests 
    pairwise_Expr= test_thisExpr_P & test_thisExpr_ES;
    pairwise_Prct= test_thisPrct_P & test_thisPrct_ES;
    
    %limiters are to force constraint on the other ES test
    if params.limitMode=="pairwise"
        pairwise_Expr=pairwise_Expr & test_pairwise_ESprct_limit;
        pairwise_Prct=pairwise_Prct & test_pairwise_ESexpr_limit;
    end
    
%     pairwise_Expr_spec=pairwise_Expr;
%     pairwise_Prct_spec=pairwise_Prct;
%     if params.specificMode=="pairwise"
%         pairwise_Expr_spec=pairwise_Expr_spec & test_thisSpecific;
%         pairwise_Prct_spec=pairwise_Prct_spec & test_thisSpecific;
%     end
    
    
    %combine all the pairwise tests. ALL by default.
    switch other.poolmethod
        case "all"
%             dominant_otherpooled(:,i)=all(pairwise_combinedES,2);
%             specific_otherpooled(:,i)=all(pairwise_combinedES_spec,2);
            self_test_exprES(:,i)=all(pairwise_Expr,2);
            self_test_prctES(:,i)=all(pairwise_Prct,2);
        case "any"
%             dominant_otherpooled(:,i)=any(pairwise_combinedES,2);
%             specific_otherpooled(:,i)=any(pairwise_combinedES_spec,2);
            self_test_exprES(:,i)=any(pairwise_Expr,2);
            self_test_prctES(:,i)=any(pairwise_Prct,2);
        case "some"
%             dominant_otherpooled(:,i)=sum(pairwise_combinedES,2)>other.some_N;
%             specific_otherpooled(:,i)=sum(pairwise_combinedES_spec,2)>other.some_N;
            self_test_exprES(:,i)=sum(pairwise_Expr,2)>=other.some_N;
            self_test_prctES(:,i)=sum(pairwise_Prct,2)>=other.some_N;
    end
    
    % combine the two ES tests. After pooling others (e.g. all exprtest OR
    % all prcttest - not a mixture of individual pairs...)
    switch params.combine_prct_expr
        case 'or'
            %marks dominant if some tests pass only by expr, others only by
            %prct such that neither all(expr tests) nor all(prct test).
%             pairwise_combinedES=pairwise_Expr|pairwise_Prct; %permits either test
%             pairwise_combinedES_spec=pairwise_Expr_spec|pairwise_Prct_spec;
            dominant_otherpooled(:,i)=self_test_exprES(:,i)|self_test_prctES(:,i);
        case 'and'
%             pairwise_combinedES=pairwise_Expr&pairwise_Prct; %strict
%             pairwise_combinedES_spec=pairwise_Expr_spec&pairwise_Prct_spec; 
            dominant_otherpooled(:,i)=self_test_exprES(:,i)&self_test_prctES(:,i);
    end
end

%apply uniform limiters and specific test (to remove stuff from "any" cases)
% if params.limitMode=="uniform"
%     dominant_otherpooled=dominant_otherpooled & test_pairwise_ESexpr_limit;
% end
% if params.specificMode=="uniform"
%     specific_otherpooled=dominant_otherpooled & all(test_thisOtherPrct,2);
% end

    
% test_prct and test_specific
dominant_otherpooled = dominant_otherpooled & test_minPrct;
specific_otherpooled = dominant_otherpooled & test_specific;

    
%pairwise self pooling
switch self.poolmethod
    case "all"
        dominant_selfpooled=all(dominant_otherpooled,2);
        specific_selfpooled=all(specific_otherpooled,2);
        exprTestPooled=all(self_test_exprES,2);
        prctTestPooled=all(self_test_prctES,2);
    case "any"
        dominant_selfpooled=any(dominant_otherpooled,2);
        specific_selfpooled=any(specific_otherpooled,2);
        exprTestPooled=any(self_test_exprES,2);
        prctTestPooled=any(self_test_prctES,2);
    case "some" 
        dominant_selfpooled=sum(dominant_otherpooled,2)>=self.some_N;
        specific_selfpooled=sum(specific_otherpooled,2)>=self.some_N;
        exprTestPooled=sum(self_test_exprES,2)>=self.some_N;
        prctTestPooled=sum(self_test_prctES,2)>=self.some_N;
end


%assemble the final selections
dominantCondition = dominant_selfpooled;
specificCondition = specific_selfpooled;

%tables for output
allgenes=table();
allgenes.id=genes.id;
allgenes.name=genes.name;


% allgenes.dominant=dominantCondition; %mainly for allgenes
allgenes.specific=specificCondition;

allgenes.min_self_prct=min(prct_self,[],2);
allgenes.max_other_prct=max(prct_other,[],2);
allgenes.min_self_expr=min(expr_self,[],2);
allgenes.max_other_expr=max(expr_other,[],2);

%effect sizes
allgenes.min_d_prct=min(pairwise_d_prct,[],2);
allgenes.min_fc_expr=min(pairwise_fc_expr,[],2);
% allgenes.min_log2fc_expr=min(pairwise_log2fc_expr,[],2);
% allgenes.max_log2fc_expr=max(pairwise_log2fc_expr,[],2);
% allgenes.min_d_prct_rel=min(pairwise_d_prct_rel,[],2);
% allgenes.max_d_prct_rel=max(pairwise_d_prct_rel,[],2);

allgenes.max_pwp=max(pairwise_Pmc,[],2);
allgenes.max_pwz=max(pairwise_Pz,[],2);

% anova, chi2, and worst pariwise p-values
allgenes.p_anova=Panova;
allgenes.min_pwp=min(pairwise_Pmc,[],2);
allgenes.p_chi2=Pchi2;
allgenes.min_pwz=min(pairwise_Pz,[],2);

allgenes.all_expr_test=exprTestPooled;
allgenes.all_prct_test=prctTestPooled;

allgenes.max_d_prct=max(pairwise_d_prct,[],2);
allgenes.max_fc_expr=max(pairwise_fc_expr,[],2);

allgenes.max_self_prct=max(prct_self,[],2);
allgenes.min_other_prct=min(prct_other,[],2);

allgenes.max_self_expr=max(expr_self,[],2);
allgenes.min_other_expr=min(expr_other,[],2);


%raw expression/prct values per type
allgenes=[allgenes,genes(:,...
    ismember(genes.Properties.VariableNames,[prctselfnames;prctothernames])...
   |ismember(genes.Properties.VariableNames,[exprselfnames;exprothernames]))];

%effect sizes
% allgenes=[allgenes,array2table(pairwise_prctES,'variablenames',strcat(prctESname,'_',fcCombs(:,1),'_',fcCombs(:,2)))];
% allgenes=[allgenes,array2table(pairwise_exprES,'variablenames',strcat(exprESname,'_',fcCombs(:,1),'_',fcCombs(:,2)))];

%pvals
% allgenes=[allgenes,array2table(pairwise_Pmc,'variablenames',strcat('pmc_',fcCombs(:,1),'_',fcCombs(:,2)))];
% allgenes=[allgenes,array2table(pairwise_Pz,'variablenames',strcat('pz_',fcCombs(:,1),'_',fcCombs(:,2)))];

if params.sortval~="none"
    [allgenes,ixs]=sortrows(allgenes,params.sortval,params.sortdir);
    dominantCondition=dominantCondition(ixs);
    specificCondition=specificCondition(ixs);
end

dominant=allgenes(dominantCondition,:);
specific=allgenes(specificCondition,:);

% if length(self.names)>1
% allgenes=allgenes(~(allgenes.max_self_expr==0)&~(allgenes.max_other_expr==0),:);
% else
% allgenes=allgenes(~(allgenes.self_expr==0)&~(allgenes.max_other_expr==0),:);
% end

end


function [selfpars,otherpars,params]=parse_input_pars(selfin,otherin, paramsin)
selfpars.names={'self'};
selfpars.minprct=0;
selfpars.maxprct=Inf;
selfpars.poolmethod='all';
selfpars.some_N=1; %as in for "any"

otherpars.names={'other'};
otherpars.minprct=0;
otherpars.maxprct=Inf;
otherpars.poolmethod='all';
otherpars.some_N=1; %as in for "any"

%these are all filtering parameters. Should be a different struct.
params.direction='any';
params.limitMode='pairwise';
params.combine_prct_expr='or';
params.specificMode='uniform'; %or pairwise
params.prctMetric='diff';
params.prctESThr=20;
params.fcExprThr=1;
params.minPrctEffect=-100;
params.minFcEffect=0;
params.pthr=0.1;
params.pthrpw=0.1;
params.pprothr=0.1;
params.pprothrpw=0.1;
params.sortval='none';
params.sortdir='descend';

%overwrite field values with any that were passed in
if nargin>0
for f=intersect(fieldnames(selfpars),fieldnames(selfin))'
    selfpars.(f{1})=selfin.(f{1});
end
end
if nargin>1
for f=intersect(fieldnames(otherpars),fieldnames(otherin))'
    otherpars.(f{1})=otherin.(f{1});
end
end
if nargin>2
for f=intersect(fieldnames(params),fieldnames(paramsin))'
    params.(f{1})=paramsin.(f{1});
end
end

% if any(~ismember(fieldnames(selfin),fieldnames(self)))
%     warning('unused fields in self-struct')
% end
% if any(~ismember(fieldnames(otherin),fieldnames(other)))
%     warning('unused fields in self-struct')
% end

end

