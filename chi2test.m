function [h,p, stats]=chi2test(counts, alpha)
%Pearson's chi-squared test
%counts is the contingency table according to two factors.

if ~exist('alpha','var')||isempty(alpha)
    alpha=0.05;
end

[nRows,nCols] = size(counts);

%marginals
colsums=sum(counts,1);
rowsums=sum(counts,2);
total=sum(counts(:));
dof=max(1,nRows-1)*max(1,nCols-1);

expected=rowsums*colsums/total;

chi2=(counts-expected).^2./expected; %use this to get p-vals for each pair -> correct for multcompare?
chi2stat=sum(chi2(:));

strength=sqrt( chi2stat /(total*min(nRows-1,nCols-1)));

p = chi2cdf(chi2stat,dof,'upper');

h = p<alpha;

stats.chi2stat=chi2stat;
stats.chi2=chi2;
stats.dof=dof;
stats.expected=expected;
stats.strength=strength;