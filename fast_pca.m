function [coeff,score]=fast_pca(X,npc,maxscaled)
% X is NxD matrix of data

%TODO: increase tolerance to increase speed?

if exist('maxscaled','var')
    X=(X-mean(X,1))./std(X,[],1);
    X(X>maxscaled)=maxscaled;
    X(X<-maxscaled)=-maxscaled;
end

% n=size(X,1);
% DOF = max(0,n-1);%expects centered data
[U, sigma, coeff] = svds(X, npc);  %irlba in MatLab
score =  U*sigma';

%also return eigenvalues of covariance matrix: si^2/(n-1) ?