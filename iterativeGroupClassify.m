function [cellID,new_ct,IDs,TF,mscores]=...
    iterativeGroupClassify(tcounts,genes,init_ct,params,scores)
%iterative classifier given initial celltype tree

%initial partition
[cellID,IDs,TF,mscores]=init_ct.classifyByScore(tcounts,genes);

if params.doImpute
    new_cellID=cellID;
    cellsToImpute=new_cellID=="Unc";
    [new_cellID_imp, newID]=imputeCellType(new_cellID,cellsToImpute,params.impute,scores);
    ix=find(cellsToImpute); 
    new_cellID_imp(ix(newID=="Amb"))="Unc"; %don't impute to "Ambiguous" class.
    numImputed=nnz(new_cellID=="Unc")-nnz(new_cellID_imp=="Unc")
    new_cellID=new_cellID_imp;
    cellID=new_cellID;
end
summary(cellID)

iter=1;
% new_ct=copy(init_ct);
new_ct=init_ct;
while iter<=params.maxIter
    
    % for each type in cellID that is not Amb/Unc, try to refine markers:
    cats=string(categories(cellID));
    cats(cats=="Amb")=[];
    cats(cats=="Unc")=[];
    for i = 1:length(cats)
        initial_markers = new_ct.subtypes(i).markers;
        %check self % and % in every other type for all genes.
        self=cellID==cats(i);
        prct_self=sum(tcounts(:,self)>genes.thr,2)./nnz(self)*100;
        keep=prct_self>params.selfPrctThr;
        
        others=setxor(cats,cats(i),'stable');
        for j=1:length(others)
            other=cellID==others(j);
            prct_other=sum(tcounts(:,other)>genes.thr,2)./nnz(other)*100;
            keep=keep & prct_other<params.otherPrctThr;
        end
        
        markers=genes.name(keep);
        if params.keep_initial_ids
            markers=unique([markers(:);initial_markers(:)]);
        end
        
        %set the new markers
        [gix,markers]=getGeneIndices(markers,genes.name);
        new_ct.subtypes(i).markers=markers;
        
        %set a new score threshold
        geneThresholds=genes.thr(gix);
        genesAbove=tcounts(gix,:)>=geneThresholds;
        markerScore=sum(genesAbove,1);
        scoreThreshold=round(geneThresholdOtsu(markerScore));
        scoreThreshold=max(init_ct.subtypes(i).threshold,scoreThreshold); %don't go below init thr
        new_ct.subtypes(i).threshold=scoreThreshold;
    end

    % classify based on new markers
    [new_cellID,IDs,TF,mscores]=new_ct.classifyByScore(tcounts,genes);
    
    if params.doImpute
        cellsToImpute=new_cellID=="Unc";
        [new_cellID_imp, newID]=imputeCellType(new_cellID,cellsToImpute,params.impute,scores);
        ix=find(cellsToImpute); 
        new_cellID_imp(ix(newID=="Amb"))="Unc"; %don't impute to "Ambiguous" class.
        numImputed=nnz(new_cellID=="Unc")-nnz(new_cellID_imp=="Unc")
        new_cellID=new_cellID_imp;
    end
    summary(new_cellID)
    
    if isequal(cellID,new_cellID)
        cellID=new_cellID;
        break
    else
        cellID=new_cellID;
    end
    
    iter=iter+1;
end