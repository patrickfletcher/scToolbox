function [classID,classImpute]=imputeCellType(classID,subset,params,coords)
% impute unclassified cells' type by looking at its kNNs' types. 
% classID should be categorical

%TODO: neighborhood voting methods?
%TODO: subset could be a name of a category in classID?

% method='mode';
% method='dWeighted';

if ~exist('subset','var')||isempty(subset)
    subset=true(size(classID));
end

if ~iscategorical(classID)
    classID=categorical(classID);
end

kNN=params.kNN;
method=params.method;

unclassified=subset;
nUnclass=nnz(unclassified);

if kNN==0 %nothing to do
    return
end

[nnix,Dix]=knnsearch(coords,coords(unclassified,:),'K',kNN+1);
% NnClass=zeros(nUnclass,kNN);
for i=1:nUnclass
    NnClass(i,:)=classID(nnix(i,2:end));
end

switch method
    case 'mode'
        classImpute=mode(NnClass,2); %assign by majority vote: mode

    case 'dist'
        %distance weighted: score_i=sum(1/distances to cells of type i), impute=max
        Dix=Dix(:,2:end);
        Wix=1./Dix;
        uNNClass=unique(NnClass(:));
        score=zeros(nUnclass,length(uNNClass));
        for i=1:nUnclass
            for j=1:length(uNNClass)
                score(i,j)=sum(Wix(i,NnClass(i,:)==uNNClass(j)));
            end
        end
        [~,maxix]=max(score,[],2);
        classImpute=uNNClass(maxix);
end

% classImpute(classImpute==-1)=0; %ambiguous is not valid class
classID(unclassified)=classImpute;

