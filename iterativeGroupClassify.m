function [cellID,new_ct,IDs,TF,mscores]=...
    iterativeGroupClassify(tcounts,genes,init_ct,params)
%iterative classifier given initial celltype tree

%initial partition
[cellID,IDs,TF,mscores]=init_ct.classifyByScore(tcounts,genes);
summary(cellID)

iter=1;
new_ct=copy(init_ct);
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
        new_ct.subtypes(i).threshold=scoreThreshold;
    end

    % classify based on new markers
    [new_cellID,IDs,TF,mscores]=new_ct.classifyByScore(tcounts,genes);
    summary(new_cellID)
    
    if isequal(cellID,new_cellID)
        cellID=new_cellID;
        break
    else
        cellID=new_cellID;
    end
    
    iter=iter+1;
end