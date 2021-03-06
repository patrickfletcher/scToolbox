function [mitoexclude, ax]=findMitoPerType(cells,group,params,figID, celltypeTree)

if ~iscategorical(group)
    group=categorical(group);
end
group=group(:);
group=removecats(group);
groupNames=categories(group)';
nTypes=length(groupNames);

mitoexclude=false(length(group),1);
for i=1:nTypes
    thisIx=group==groupNames(i);
    thisFracMT=cells.fracMT(thisIx);
    
    %mean+k*std
    centralFrac(i)=mean(thisFracMT);
%     centralFrac(i)=median(thisFracMT);
    spreadMT(i)=std(thisFracMT);
%     spreadMT(i)=mad(thisFracMT,1);

    high_thr(i) = min(centralFrac(i)+params.n_dev*spreadMT(i), 1)
    low_thr(i) = max(centralFrac(i)-params.n_dev*spreadMT(i), 0)

    switch params.tailoption
        case "high"
            mitoexclude=mitoexclude | (thisIx & cells.fracMT>high_thr(i));
        case "low"
            mitoexclude=mitoexclude | (thisIx & cells.fracMT<low_thr(i));
        case "both"
            mitoexclude=mitoexclude | (thisIx & cells.fracMT>high_thr(i));
            mitoexclude=mitoexclude | (thisIx & cells.fracMT<low_thr(i));
        otherwise 
            error('unknown tail option')
    end
end

if exist('figID','var') && ~isempty(figID)
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
    
%     x = cells.genesPerCell; xlab='Genes per GEM';
    x = cells.molecPerCell; xlab='UMI per GEM';
    
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
    
    for i=1:nsub
        ax(i)=tight_subplot(nr,nc,i,gap,marg_h,marg_w);
%         ax(i).XScale='log';

        thisIx=group==groupNames{i};

        plot(x(thisIx),cells.fracMT(thisIx),'color',colors(i,:),'linestyle','none','marker',markerlist(i))
        line(x(thisIx & mitoexclude),cells.fracMT(thisIx & mitoexclude),'color',colors(i,:)/4,'linestyle','none','marker',markerlist(i))  
        line(xlim(),centralFrac(i)*[1,1],'color',colors(i,:)/2,'linestyle','--')
        switch params.tailoption
            case "high"
                line(xlim(),high_thr(i)*[1,1],'color',colors(i,:)/4)
            case "low"
                line(xlim(),low_thr(i)*[1,1],'color',colors(i,:)/4)
            case "both"
                line(xlim(),high_thr(i)*[1,1],'color',colors(i,:)/4)
                line(xlim(),low_thr(i)*[1,1],'color',colors(i,:)/4)
        end
%         axis(ax(i),'tight')
        MING=min(x(thisIx));
        MAXG=max(x(thisIx));
        MINF=min(cells.fracMT(thisIx));
        MAXF=max(cells.fracMT(thisIx));
        ax(i).XLim=[MING*0.9,MAXG*1.1];
%         ax(i).YLim(2)=ceil(100*MAXF)/100;
        ax(i).YLim(2)=MAXF+0.03;
        ax(i).XTick=[MING,MAXG];
%         ax(i).YTick=0:0.1:0.4;
%         ax(i).YTick=[0,floor(10*MAXF)/10];
        

        if i==1
            ylabel('f_{MT}')
        end
        if i==nTypes
            xlabel(xlab)
        end
%         if i==(nr-1)*nc+1
%             xlabel('Genes per gem')
%             ylabel('Fraction mito')
%         end
        
        if exist('names','var')
            title(names(i),'FontSize',8)
        end
        box on
        ax(i).TickLength=[.05,.05];
        ax(i).FontSize=8;
    end
%     
%     ax(1).XTick=[0,500,1000];
%     for i=2:9
%     ax(i).XTick=0:3000:10000;
%     end
%     ax(4).XTick=0:4000:10000;
    
    figID.PaperUnits='centimeters';
    figID.PaperPosition=[0,0,figID.Position(3:4)];
    figID.PaperSize=figID.Position(3:4);
end