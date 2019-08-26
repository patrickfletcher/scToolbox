function X=winsorize(X,ptiles)
% clamp outliers in each row (variable) to percentiles of the distribution

PT=prctile(X,ptiles,2);

Xlo=X<PT(:,1);
Xhi=X>PT(:,2);

for i=1:size(X,1)
    X(Xlo(i,:))=PT(i,1);
    X(Xhi(i,:))=PT(i,2);
end