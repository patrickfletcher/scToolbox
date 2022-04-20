function [V, score, mu, sig]=fast_pca(X,npc,maxscaled)
% X is NxD matrix of data

%covariance = X'*X/(n-1) <=> X is centered.

% https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca

%TODO: support returning eigenvalues (PC variances)  si^2/(n-1)
%TODO: standardized scores = sqrt(n-1)*U
%TODO: clarify loadings vs V
%TODO: add interface to numerical options of eigs (tolerance, etc)

mu=zeros(size(X,1),1);
sig=ones(size(X,1),1);
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
[U,S,V] = svdsecon(X,npc);
score =  U*S';
