function [p,PZ,stats,combs]=pairwiseProportionTest(X,group)
% Perform Pearson's chi square test of independence for counts partitioned
% by group according to a binary factor, then test difference between all
% pairs of groups.

%percent expressing can be relative to threshold

%TODO: pass in X=tcounts-thr, check >=0 (remove need for genes arg)
%TODO: chi2 table vs vector?? different test statistic for the same conceptual test??

correctionmethod='fdr';

nGenes=size(X,1);

if ~iscategorical(group)
    group=categorical(group);
end
groupNames=categories(group);
groupCounts=countcats(group);
totalCells=sum(groupCounts);
nGroups=length(groupNames);

%null distribution: expectedCounts - evenly spread total counts across celltypes
observedCounts=zeros(nGenes,nGroups);
for i=1:nGroups
    g=group==groupNames{i};
    observedCounts(:,i)=sum(X(:,g)>0,2);   %or =prct/100*groupCounts
end
expectedCounts=sum(observedCounts,2)./totalCells .*groupCounts'; %probability * N_celltype
dof=nGroups-1;

% full table version - also use the <thr counts (makes things more
% significant....)
observedCounts=[observedCounts,groupCounts'-observedCounts];
expectedCounts=[expectedCounts,groupCounts'-expectedCounts];
% dof=max(1,nRows-1)*max(1,nCols-1); %nRows=nGroups; nCols=2 (>thr or <=thr)

%chi2 statistic
chi2=(observedCounts-expectedCounts).^2./expectedCounts;
chi2stat=sum(chi2,2);
p = chi2cdf(chi2stat,dof,'upper');

%if not enough expected counts, the chi2 test isn't valid...
nExpectedLT5=sum(expectedCounts<5,2);
valid5=nExpectedLT5<0.2*nGroups;
nExpectedLT1=sum(expectedCounts<1,2);
valid1=nExpectedLT1==0;

p(~valid5)=nan; 

%now do z-test for difference in proportion, pairwise
pooledVariance=false;
combs=combnk(groupNames,2);
nPairs=size(combs,1);
PZ=nan(nGenes,nPairs);
for i=1:nPairs
    selfIdx=strcmp(groupNames,combs{i,1});
    otherIdx=strcmp(groupNames,combs{i,2});
    [~,pz] = ztestp(observedCounts(:,selfIdx),groupCounts(selfIdx),observedCounts(:,otherIdx),groupCounts(otherIdx),pooledVariance);
    
    switch lower(correctionmethod)
        case 'bonferroni'
            %sequential bonferroni-holm
            PZ(:,i)=bonf_holm(pz);
        case {'fdr','bh'}
            %Benjamini-Hochberg FDR
            [~, ~, ~, PZ(:,i)]=fdr_bh(pz); %method?
        otherwise
            PZ=pz;
    end
end



PZ(PZ>1)=1;

stats.chi2stat=chi2stat;
stats.dof=dof;
stats.observedCounts=observedCounts;
stats.expectedCounts=expectedCounts;
stats.valid5=valid5;
stats.valid1=valid1;