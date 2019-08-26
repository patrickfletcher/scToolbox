function [p,h]=findSignificantPCs(data,nReps,alpha,maxPCs)
%run a permutation test on the data to determine the number of PCs to use
%in data reduction. observations in rows, variables in columns

[n,m]=size(data);
% maxPCs=min(n,m);
% maxPCs=m;


%observed svds
% tic
svsObs=svd(data);
% svsObs=svds(data,maxPCs);
% toc
% lamObs=svsObs.^2./(n-1); %variances of each component

% permutation test: shuffle labels of observations within each variable independently
dataPerm=zeros(n,m);
svs=zeros(nReps,length(svsObs));
tenth=round(nReps/10);
for it=1:nReps
    for j=1:m
        dataPerm(:,j)=data(randperm(n),j);
    end
    
%     tic
    svs(it,:)=svd(dataPerm);
%     svs(it,:)=svd(dataPerm,'econ');
%     svs(it,:)=svds(dataPerm,1);
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