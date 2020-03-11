function [npc,p,eObs,max_e_null]=findSignificantPCs(data,nReps,alpha)
%run a permutation test on the data to determine the number of PCs to use
%in data reduction. observations in rows, variables in columns

rng('shuffle','simdTwister') %for speedup?

[m,n]=size(data);

%observed svds
tic
eObs = eig(data'*data,'nobalance');
eObs = sort(abs(eObs),'descend');
% svsObs=sqrt(eObs); %actually no need to sqrt...
% svsObs=svd(data);
% svsObs=svds(data,maxPCs);
toc
% lamObs=svsObs.^2./(n-1); %variances of each component

% permutation test: shuffle labels of observations within each variable independently
dataPerm=zeros(m,n);
e_null=zeros(length(eObs),nReps);
% svs=zeros(nReps,length(svsObs));
% svs=zeros(nReps,maxPCs);
tenth=round(nReps/10);
for it=1:nReps
%     tic
    for j=1:n
        dataPerm(:,j)=data(randperm(m,m),j);
    end
%     toc %10X less time than svd below
    
%     tic
    d = eig(dataPerm'*dataPerm,'nobalance');
    e_null(:,it) = sort(abs(d),'descend');
%     d = sort(abs(d),'descend');
%     svs(it,:)=sqrt(d); %actually no need to sqrt...
%     toc

    if mod(it,tenth)==0
        fprintf('.')
    end
end
% lam=svs.^2./(n-1);

%global threshold
%Moffitt (science) threshold: mean of distribution of null max eig -  svsObs>mean(svs(:,1))

%rank observed eigs into null distribution of largest eigs. From above,
%largest eigenvalue in each trial is in row 1 of e_null.
max_e_null=sort(e_null(1,:),'descend');
ranks=zeros(size(eObs));
for i=1:length(eObs)
    firstLower=find(max_e_null<eObs(i),1,'first');
    if ~isempty(firstLower)
        ranks(i)=firstLower;
    else
        ranks(i)=nReps;
    end
end
p=ranks/nReps;

%find the number of components based on alpha
h=p<alpha;
nPCs=find(~h,1,'first')-1;
    
% fprintf('\n')

% %rank observed SVs into distribution of null largest svs
% svs1Sorted=sort(svs(:,1),'descend');
% ranks=zeros(size(svsObs));
% for i=1:length(svsObs)
%     firstLower=find(svs1Sorted<svsObs(i),1,'first');
%     if ~isempty(firstLower)
%         ranks(i)=firstLower;
%     else
%         ranks(i)=nReps;
%     end
% end
% p=ranks/nReps;

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