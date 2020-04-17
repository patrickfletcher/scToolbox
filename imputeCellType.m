function [classID,classImpute,to_impute]=imputeCellType(classID,subset,params,coords)
% impute unclassified cells' type by looking at its kNNs' types. 
% classID should be categorical

%TODO: allow for passing in pre-computed NNs
%TODO: neighborhood voting methods? probability based interpretation?
%TODO: subset could be a name of a category in classID?

% method='mode';
% method='dWeighted';


if ~iscategorical(classID)
    classID=categorical(classID);
end

kNN=params.kNN;
method=params.method;

to_impute=true(size(classID));
if ~isempty(subset)
    if isstring(subset)||ischar(subset)||iscellstr(subset)
        to_impute=classID==subset;
    elseif size(subset)==size(classID)
        to_impute=subset;
    end
end
n_to_impute=nnz(to_impute);

if kNN==0 %nothing to do
    return
end

[nnix,Dix]=knnsearch(coords,coords(to_impute,:),'K',kNN+1);

for i=1:n_to_impute
    NnClass(i,:)=classID(nnix(i,:));
end

switch method
    case 'mode'
        classImpute=mode(NnClass,2); %assign by majority vote: mode

    case 'nearest'
        %nearest neighbor that is not also unc
        classImpute=classID(to_impute)';
        for i=1:n_to_impute
            myID=NnClass(i,1);
            thisnix=find(NnClass(i,:)~=myID&NnClass(i,:)~="Unc",1,'first'); %what if this fails? returns []
            if ~isempty(thisnix)
                nix(i,1)=thisnix;
                d(i,1)=Dix(i,nix(i));
                classImpute(i,1)=NnClass(i,nix(i));
            end %otherwise, stays unc.
        end
        
    case 'dist'
        %distance weighted: score_i=sum(1/distances to cells of type i), impute=max
        Dix=Dix(:,2:end);
        Wix=1./Dix;
        uNNClass=unique(NnClass(:));
        score=zeros(n_to_impute,length(uNNClass));
        for i=1:n_to_impute
            for j=1:length(uNNClass)
                score(i,j)=sum(Wix(i,NnClass(i,:)==uNNClass(j)));
            end
        end
        [~,maxix]=max(score,[],2);
        classImpute=uNNClass(maxix);
end

% classImpute(classImpute==-1)=0; %ambiguous is not valid class
classID(to_impute)=classImpute;

