function [coeff,score, mu, sig]=fast_pca(X,npc,maxscaled)
% X is NxD matrix of data

%TODO: increase tolerance to increase speed?

mu=zeros(size(X,1),1);
sig=zeros(size(X,1),1);
if exist('maxscaled','var')
    mu=mean(X,1);
    sig=std(X,[],1);
    X=(X-mu)./sig;
    X(X>maxscaled)=maxscaled;
    X(X<-maxscaled)=-maxscaled;
end

% n=size(X,1);
% DOF = max(0,n-1);%expects centered data
% [U, sigma, coeff] = svds(X, npc);  %irlba in MatLab
[U,S,coeff] = svdsecon(X,npc);
score =  U*S';

%also return eigenvalues of covariance matrix: si^2/(n-1) ?