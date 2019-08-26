function T=boxcox(X,lam)
%rows are variables, colums observations

H=geomean(X,2);
if lam==0
    T=H.*log(X);
else
    T=(X.^lam - 1)./(lam*H.^(lam-1));
end