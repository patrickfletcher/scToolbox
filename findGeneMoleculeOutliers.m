function [outliers,ax, ptilesGenes,ptilesMolec]=findGeneMoleculeOutliers(cells, group, params, figID, celltypeTree)

% pct specifies total percent of cells to exclude; if both tails, get half
% lowest half highest.

if ~iscategorical(group)
    group=categorical(group);
end
group=group(:);
group=removecats(group);
groupNames=categories(group)';
nTypes=length(groupNames);

outliers=false(size(cells.genesPerCell));
for i=1:nTypes
    thisIx=group==groupNames(i);
    thisGenes=cells.genesPerCell(thisIx);
    thisMolec=cells.molecPerCell(thisIx);
    
    switch params.method
        case "pctile"
            if strcmpi(params.tailoption,'both')
                pct=pct/2;
            end
            ptilesGenes(i,:)=prctile(thisGenes,[pct,100-pct]);
            gene_low=ptilesGenes(i,1);
            gene_high=ptilesGenes(i,2);
            ptilesMolec(i,:)=prctile(thisMolec,[pct,100-pct]);
            molec_low=ptilesMolec(i,1);
            molec_high=ptilesMolec(i,2);
        case "values"
            gene_low=params.min_genes;
            gene_high=params.max_genes;
            molec_low=params.min_molec;
            molec_high=params.max_molec;
    end
    
    
    
    switch lower(params.tailoption)
        case 'high'
            outlierGenes{i}=cells.genesPerCell>gene_high;
            outlierMolecs{i}=cells.molecPerCell>molec_high;
        case 'low'
            outlierGenes{i}=cells.genesPerCell<gene_low;
            outlierMolecs{i}=cells.molecPerCell<molec_low;
        case 'both'
            outlierGenes{i}=cells.genesPerCell<gene_low|cells.genesPerCell>gene_high;
            outlierMolecs{i}=cells.molecPerCell<molec_low|cells.molecPerCell>molec_high;
        otherwise 
            error('unknown tail option')
    end
    
    outliers=outliers|(thisIx & outlierGenes{i});
    outliers=outliers|(thisIx & outlierMolecs{i});
end

%TODO: move plotting to separate function?

if exist('figID','var') & ~isempty(figID)
    if exist('celltypeTree','var')
        colors=celltypeTree.Colors(groupNames);
        names=celltypeTree.Names(groupNames);
    else
        if nTypes==1
            colors=[0.5,0.5,0.5];
        else
            colors = cbrewer('qual','Set1',min(3,nTypes));
        end
    end

    markerlist=repmat('.',1,nTypes);
%     markerlist='+o*.xsdph^v><';
    
    figID=figure(figID);clf
    figID.Units='centimeters';
    figID.Position(3)=10;
    figID.Position(4)=7.5;

    gap=[0.1,0.1];
    marg_h=[0.125,0.075];
    marg_w=[0.11,0.04];
    
    nsub=nTypes;

    p=numSubplots(nsub);
    nr=p(1); nc=p(2);
    
    cells.molecPerCell=cells.molecPerCell/10000;
    for i=1:nsub
%         ax(i)=subplot(nr,nc,i);
        ax(i)=tight_subplot(nr,nc,i,gap,marg_h,marg_w);
    
        thisIx=group==groupNames(i);
        line(cells.genesPerCell(thisIx),cells.molecPerCell(thisIx),'color',colors(i,:),'marker',markerlist(i),'linestyle','none')

        
%         line(medGenes(i)*[1,1],ylim(),'color',colors(i,:),'linestyle','--')
%         line(xlim(),medMolec(i)*[1,1],'color',colors(i,:),'linestyle','--')

        line(cells.genesPerCell(thisIx & outlierGenes{i}),cells.molecPerCell(thisIx & outlierGenes{i}),'color',colors(i,:)/4,'marker',markerlist(i),'linestyle','none')
        line(cells.genesPerCell(thisIx & outlierMolecs{i}),cells.molecPerCell(thisIx & outlierMolecs{i}),'color',colors(i,:)/4,'marker',markerlist(i),'linestyle','none')

%         XLIM=xlim(); YLIM=ylim();
        MING=min(cells.genesPerCell(thisIx));
        MAXG=max(cells.genesPerCell(thisIx));
        MINM=min(cells.molecPerCell(thisIx));
        MAXM=max(cells.molecPerCell(thisIx));
        ax(i).XLim=[MING*0.9,MAXG*1.1];
        ax(i).YLim=[MINM*0.9,MAXM*1.1];
        ax(i).XTick=[MING,MAXG];
        ax(i).YTick=[ceil(MINM*100),floor(MAXM*100)]./100;

%         ax(i).XScale='log';
%         ax(i).YScale='log';
        
        if i==1
            ylabel('UMI per GEM (x10^4)')
            ax(i).YTick=[ceil(MINM*100),floor(MAXM/2*100),floor(MAXM*100)]./100;
        end
        if i==nTypes
            xlabel('Genes per GEM')
        end
%         if i==(nr-1)*nc+1
%             xlabel('Genes per GEM')
%             ylabel('UMI per GEM')
%         end
        if exist('names','var')
            title(names(i),'FontSize',8)
        end
        box on
        ax(i).TickLength=[.05,.05];
        ax(i).FontSize=8;
    end
    
    figID.PaperUnits='centimeters';
    figID.PaperPosition=[0,0,figID.Position(3:4)];
    figID.PaperSize=figID.Position(3:4);
end
