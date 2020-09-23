function [Y,cellPerm]=groupCountMatrix(X,ident,orderCellsOption,PCscores)

minK=10;

if ~iscategorical(ident)
    ident=categorical(ident);
end
ident=ident(:)'; %force row vector
groupNames=categories(ident);
groupCounts=countcats(ident);

doClustOnPCscores=false;
if exist('PCscores','var')
    doClustOnPCscores=true;
end


%for subclustering:
distanceMetric='Euclidean';
linkageType='Ward';

Y=[];
cellPerm=[];
for i=1:length(groupNames)
    thisGroup=find(ident==groupNames(i));
    if groupCounts(i)>0
        
        Ysub=X(:,thisGroup);
        
        subOrder=1:groupCounts(i);
        
        %reorder cells within the group for better visual effect
        if groupCounts(i)>1

            switch orderCellsOption
                case 'none'
                    ixs=1:size(Ysub,2);
                    
                case 'nnz'
                    NNZ=sum(Ysub>0,1);
                    [~,ixs]=sort(NNZ,'descend');
                    
                case 'mean'
                    MEAN=mean(Ysub,1);
                    [~,ixs]=sort(MEAN,'descend');
        
                case 'optim'
%                     K=min(groupCounts(i),minK);
                    if doClustOnPCscores
                        xD=pdist(PCscores(thisGroup,:),distanceMetric); 
                    else
                        xD=pdist(Ysub',distanceMetric); 
                    end
                    Z=linkage(xD,linkageType);
%                     T=cluster(Z,'maxClust',K);
    
                    %slow when nObs is large
                    ixs=optimalleaforder(Z,xD);
                    % ixs=optimalleaforder(Z,xD,'Transformation','linear');
%                     ixs=optimalleaforder(Z,xD,'Transformation','inverse');

                case 'clust'
                    K=min(groupCounts(i),minK);
                    if doClustOnPCscores
                        xD=pdist(PCscores(thisGroup,:),distanceMetric); 
                    else
                        xD=pdist(Ysub',distanceMetric); 
                    end
                    Z=linkage(xD,linkageType);
                    T=cluster(Z,'maxClust',K);
                    
                    [~,ixs]=sort(T);

                otherwise % 'none'
                    ixs=1:size(Ysub,2);
                    
            end
            subOrder=subOrder(ixs);
        end
        
        Ysub=Ysub(:,subOrder);
        
        Y=[Y,Ysub];
        cellPerm=[cellPerm,thisGroup(subOrder')];
    end
end
