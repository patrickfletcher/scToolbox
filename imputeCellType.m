function [output_class, class_imputed, to_impute]=imputeCellType(input_class,subset,params,coords)
% impute unclassified cells' type by looking at its kNNs' types. 
% classID should be categorical

%TODO: allow for passing in pre-computed NNs
%TODO: neighborhood voting methods? probability based interpretation?
%TODO: subset could be a name of a category in classID?

if ~iscategorical(input_class)
    input_class=categorical(input_class);
end

if ~isnumeric(params.n_neighbors) && params.n_neighbors=="sqrtN"
    params.n_neighbors=floor(sqrt(length(input_class)));
end

if ~isfield(params,'metric')
    params.metric='euclidean';
end

result=params;

to_impute=true(size(input_class));
if ~isempty(subset)
    if isstring(subset)||ischar(subset)||iscellstr(subset)
        to_impute=input_class==subset;
    elseif size(subset)==size(input_class)
        to_impute=subset;
    end
end
n_to_impute=nnz(to_impute);

[nnix,Dix]=knnsearch(coords,coords(to_impute,:),'K',params.n_neighbors+1,'Distance',params.metric);
nnix(:,1)=[]; %remove self
Dix(:,1)=[];

for i=1:n_to_impute
    NnClass(i,:)=input_class(nnix(i,:));
end

switch params.method
    case 'mode'
        class_imputed=mode(NnClass,2); %assign by majority vote: mode

    case 'mode_notme'
        %nearest neighbor that is not also unc
        class_imputed=input_class(to_impute)';
        for i=1:n_to_impute
            myID=NnClass(i,1);
            thisnns=NnClass(i,NnClass(i,:)~=myID & NnClass(i,:)~="Unc");
            if ~isempty(thisnns)
                class_imputed(i,1)=mode(thisnns);
            end %otherwise, stays unc.
        end
        
    case 'nearest'
        %nearest neighbor that is not also unc
        class_imputed=input_class(to_impute)';
        for i=1:n_to_impute
            myID=NnClass(i,1);
            thisnix=find(NnClass(i,:)~=myID&NnClass(i,:)~="Unc",1,'first'); %what if this fails? returns []
            if ~isempty(thisnix)
                class_imputed(i,1)=NnClass(i,thisnix);
%                 nix(i,1)=thisnix;
%                 d(i,1)=Dix(i,nix(i));
%                 classImpute(i,1)=NnClass(i,nix(i));
            end %otherwise, stays unc.
        end
        
    case 'dist'
        %distance weighted: score_i=sum(1/distances to cells of type i), impute=max
        Wix=1./Dix;
        uNNClass=unique(NnClass(:));
        score=zeros(n_to_impute,length(uNNClass));
        for i=1:n_to_impute
            for j=1:length(uNNClass)
                score(i,j)=sum(Wix(i,NnClass(i,:)==uNNClass(j)));
            end
        end
        [~,maxix]=max(score,[],2);
        class_imputed=uNNClass(maxix);
end

% classImpute(classImpute==-1)=0; %ambiguous is not valid class
output_class=input_class;
output_class(to_impute)=class_imputed;

% result.output_class=output_class;
% result.class_imputed=class_imputed;
% result.to_impute=to_impute;
