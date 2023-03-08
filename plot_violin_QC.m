function [fh,ax]=plot_violin_QC(QCdata, TF, sampleID, ctnames, cols, figpos)
arguments
    QCdata
    TF 
    sampleID
    ctnames 
    cols 
    figpos=[]
end

samples=categories(sampleID);

thrw=0.33;
datasz=1;
%'Units','inches','Position',[2,2,4,2]
gap=0.01;
marg_w=[0.15,0.1];
marg_h=[0.1,0.05];

% cols=raw.cols;
% cols(end-1,:)=[]; %remove Amb
% cols(1,:)=[]; %remove E

% ylabs

% TF=raw.TF_imp;
% TF(end+1,:)=~any(TF,1);

% qcctnames=ctnames;
% qcctnames(end+1)={'Unc'};  %Amb are hidden in the other types
% qcctnames(1)=[]; %remove E
% TF(1,:)=[];


QCdataFields=fieldnames(QCdata(1));
QCdataFields(QCdataFields=="samplename")=[];
ylabs = QCdataFields;

nplots=length(QCdata);
for pix=1:nplots
    
    fh(pix)=figure();clf
    if ~isempty(figpos)
    %     fh(pix).Units='inches';
        fh(pix).Position=figpos;
    %     fh(pix).PaperPositionMode='auto'; 
    end
    nPanels=length(QCdataFields);
    for i=1:nPanels
        ax(i)=tight_subplot(nPanels,1,i,gap,marg_h,marg_w);
        thisQC=QCdata(pix).(QCdataFields{i});
        thisvals=thisQC.vals;
        thislow=thisQC.lowthr;
        thishi=thisQC.hithr;
        if thisQC.params.logval
            thisvals=log10(thisvals);
%             thislow=log10(thislow);
%             thishi=log10(thishi);
        end
        thissamplecells=sampleID==samples{pix};
        for j=1:size(TF,1)
            thisTF=TF(j,thissamplecells);
            valsub=thisvals(thisTF);
            
            %don't plot thresholds for low num cts
            if length(valsub)<30
                continue
            end
            
            hvn(j)=Violin(valsub, j,'ViolinColor',cols(j,:));
            hvn(j).MedianPlot.Marker='+';
            hvn(j).ScatterPlot.SizeData=datasz;
%             hvn(j).ShowData=false;
            hvn(j).ViolinAlpha=1;
%             if j==length(qcctnames)
            hvn(j).MedianPlot.Visible=false;
            hvn(j).BoxPlot.Visible=false;
%             end
            
            if thishi(j)>=min(thisvals)
                line(ax(i),j+thrw*[-1,1],thislow(j)*[1,1],'color','k')
            end
            if thishi(j)<=max(thisvals)
                line(ax(i),j+thrw*[-1,1],thishi(j)*[1,1],'color','k')
            end
        end
        
        if i<nPanels
            ax(i).XTick=[];
        else
            ax(i).XTick=1:length(ctnames);
        end
        ax(i).XTickLabel=ctnames;
        ax(i).XLim=[0.5,length(ctnames)+0.5];
%         if mod(i,2)==0
%             ax(i).YAxisLocation='right';
%         end
        ylabel(ax(i),ylabs{i})
        
        if thisQC.params.logval
            ax(i).YTick=0:4;
            ax(i).YTickLabels={'10^0','10^1','10^2','10^3','10^4'};
        else
            thishi(thishi==thisQC.params.hithr_clipval)=[];
            ax(i).YLim=[0,min(1.5*max(thishi),1)];
        end
    end
end