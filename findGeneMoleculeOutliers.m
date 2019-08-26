function [outliers,ptilesGenes,ptilesMolec]=findGeneMoleculeOutliers(A, group, pct, tailOption, figID, celltypeTree)

doPlot=exist('figID','var') & ~isempty(figID);

genesPerGem=sum(A>0,1);  %median of this gives value in web summary
molecPerGem=sum(A,1);

% pct specifies total percent of cells to exclude; if both tails, get half
% lowest half highest.
if strcmpi(tailOption,'both')
    pct=pct/2;
end

if ~iscategorical(group)
    group=categorical(group);
end
group=removecats(group);
groupNames=categories(group)';
nTypes=length(groupNames);

outliers=false(size(genesPerGem));
for i=1:nTypes
    thisIx=group==groupNames(i);
    thisGenes=genesPerGem(thisIx);
    thisMolec=molecPerGem(thisIx);
    
    medGenes(i)=median(thisGenes);
    medMolec(i)=median(thisMolec);
    ptilesGenes(i,:)=prctile(thisGenes,[pct,100-pct]);
    ptilesMolec(i,:)=prctile(thisMolec,[pct,100-pct]);
    
    switch lower(tailOption)
        case 'left'
            outlierGenes{i}=genesPerGem<ptilesGenes(i,1);
            outlierMolecs{i}=molecPerGem<ptilesMolec(i,1);
        case 'right'
            outlierGenes{i}=genesPerGem>ptilesGenes(i,2);
            outlierMolecs{i}=molecPerGem>ptilesMolec(i,2);
        case 'both'
            outlierGenes{i}=genesPerGem<ptilesGenes(i,1)|genesPerGem>ptilesGenes(i,2);
            outlierMolecs{i}=molecPerGem<ptilesMolec(i,1)|molecPerGem>ptilesMolec(i,2);
        otherwise 
            error('unknown tail option')
    end
    
    outliers=outliers|(thisIx & outlierGenes{i});
    outliers=outliers|(thisIx & outlierMolecs{i});
end

%TODO: move plotting to separate function?

if doPlot
    colors=celltypeTree.Colors(groupNames);
    names=celltypeTree.Names(groupNames);

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
    
    molecPerGem=molecPerGem/10000;
    for i=1:nsub
%         ax(i)=subplot(nr,nc,i);
        ax(i)=tight_subplot(nr,nc,i,gap,marg_h,marg_w);
    
        thisIx=group==groupNames(i);
        line(genesPerGem(thisIx),molecPerGem(thisIx),'color',colors(i,:),'marker',markerlist(i),'linestyle','none')

        
%         line(medGenes(i)*[1,1],ylim(),'color',colors(i,:),'linestyle','--')
%         line(xlim(),medMolec(i)*[1,1],'color',colors(i,:),'linestyle','--')

        line(genesPerGem(thisIx & outlierGenes{i}),molecPerGem(thisIx & outlierGenes{i}),'color',colors(i,:)/4,'marker',markerlist(i),'linestyle','none')
        line(genesPerGem(thisIx & outlierMolecs{i}),molecPerGem(thisIx & outlierMolecs{i}),'color',colors(i,:)/4,'marker',markerlist(i),'linestyle','none')

%         XLIM=xlim(); YLIM=ylim();
        MING=min(genesPerGem(thisIx));
        MAXG=max(genesPerGem(thisIx));
        MINM=min(molecPerGem(thisIx));
        MAXM=max(molecPerGem(thisIx));
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
        title(names(i),'FontSize',8)
        box on
        ax(i).TickLength=[.05,.05];
        ax(i).FontSize=8;
    end
    
    figID.PaperUnits='centimeters';
    figID.PaperPosition=[0,0,figID.Position(3:4)];
    figID.PaperSize=figID.Position(3:4);
end
