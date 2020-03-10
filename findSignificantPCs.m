function [p,h]=findSignificantPCs(data,nReps,alpha,maxPCs)
%run a permutation test on the data to determine the number of PCs to use
%in data reduction. observations in rows, variables in columns

rng('shuffle','simdTwister') %for speedup?

[n,m]=size(data);
% maxPCs=min(n,m);
% maxPCs=m;

% maxPCs=100;

%observed svds
% tic
e = eig(data'*data);
svsObs=sqrt(sort(abs(e),'descend'));
% svsObs=svd(data);
% svsObs=svds(data,maxPCs);
% toc
% lamObs=svsObs.^2./(n-1); %variances of each component

% permutation test: shuffle labels of observations within each variable independently
dataPerm=zeros(n,m);
svs=zeros(nReps,length(svsObs));
% svs=zeros(nReps,maxPCs);
tenth=round(nReps/10);
for it=1:nReps
%     tic
    for j=1:m
        dataPerm(:,j)=data(randperm(n,n),j);
    end
%     toc %10X less time than svd below
    
%     tic
    d = eig(dataPerm'*dataPerm);
    svs(it,:) = sqrt(sort(abs(d),'descend')); 
%     svs(it,:)=svd(dataPerm);
%     svs(it,:)=svdecon(dataPerm);
%     svs(it,:)=svd(dataPerm,'econ');
%     svs(it,:)=svd(dataPerm,0);
%     svs(it,:)=svds(dataPerm,maxPCs);
%     toc

    if mod(it,tenth)==0
        fprintf('.')
    end
end
% lam=svs.^2./(n-1);

%global threshold
%Moffitt (science) threshold: mean of distribution of null max eig -  svsObs>mean(svs(:,1))

% sv>max(svs(:,1))
%rank observed SVs into distribution of null largest svs
svs1Sorted=sort(svs(:,1),'descend');
ranks=zeros(size(svsObs));
for i=1:length(svsObs)
    firstLower=find(svs1Sorted<svsObs(i),1,'first');
    if ~isempty(firstLower)
        ranks(i)=firstLower;
    else
        ranks(i)=nReps;
    end
end
p=ranks/nReps;

% %rank observed SVs into null distributions
% ranks=zeros(size(svsObs));
% for i=1:length(svsObs)
%     svsiSorted=sort(svs(:,i),'descend');
%     firstLower=find(svsiSorted<svsObs(i),1,'first');
%     if ~isempty(firstLower)
%         ranks(i)=firstLower;
%     else
%         ranks(i)=nReps;
%     end
% end
% p=ranks/nReps;

%find the number of components based on alpha
h=p<alpha;

fprintf('\n')