function [h,p,stats] = ztestp(n1,N1,n2,N2,pooled,alpha,tails)

%sign of Z tells if p1>p2 or vice versa. 

if ~exist('pooled','var')||isempty(pooled)
    pooled=true;
end
if ~exist('alpha','var')||isempty(alpha)
    alpha=0.05;
end
if ~exist('tails','var')||isempty(tails)
    tails='both';
end

p1hat=n1./N1;
p2hat=n2./N2;

if pooled
    phat=(n1+n2)./(N1+N2);
    Z=(p1hat-p2hat)./sqrt(phat.*(1-phat).*(1/N1+1/N2)); %pooled
else
    Z=(p1hat-p2hat)./sqrt( p1hat.*(1-p1hat)./N1 + p2hat.*(1-p2hat)./N2); %unpooled - different variances
end

%handle upper/lower/two-tail
switch tails
    case 'upper'
        p=normcdf(Z,'upper');
        h=p<alpha;
    case 'lower'
        p=normcdf(Z);
        h=p<alpha;
    case 'both'
        p=2*normcdf(abs(Z),'upper');
        h=p<alpha;
%         p=normcdf(abs(Z),'upper');
%         h=p<alpha/2;
end

stats.Z=Z;
stats.estimate=p1hat-p2hat;