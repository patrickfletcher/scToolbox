function hbc=clustered_bubblechart(termtab, groupnames, sizevar, colvar, groupvar, options)
arguments
    termtab
    groupnames
    sizevar='p_value'
    colvar='fold_enrichment'
    groupvar=[];
    options.ysortvals='min_p_value'
    options.ysortval_dir='ascend'
    options.max_termlength=40
    options.ylabgroup=[] %termtab variable
    options.ylabcol=[]
    options.show_source=true;
    options.show_ygroup_label=false
    options.axpos=[0.125,0.05,0.25,0.85]
    options.sizerange=[4,12]
    options.sizelims=[eps,1];
    options.colrange=[]
    options.group_order=[]
    options.log_color=true
    options.log_color_mode='two'
    options.cb_ticks=3
    options.log_size=true
    options.log_size_mode='-ten'
    options.sigfigs=2
end
%only plotting. do selection outside.

sortcols={};
if isempty(groupvar)
    groupvar=ones(height(termtab),1);
    nPerGroup=[];
else
    nPerGroup=countcats(categorical(termtab.(groupvar)));
    sortcols={groupvar};
end

if ~isempty(options.ylabgroup)
    sortcols(end+1)={options.ylabgroup};
end

sortcols(end+1)={options.ysortvals};
sordir=options.ysortval_dir;
%TODO: sorting is disabled
% termtab=sortrows(termtab,sortcols,sordir);

sizlab=strrep(sizevar,'_',' ');
SIZ=termtab{:,startsWith(termtab.Properties.VariableNames, sizevar)};

%size limits?
SIZ(SIZ<options.sizelims(1))=options.sizelims(1);
SIZ(SIZ>options.sizelims(2))=options.sizelims(2);

collab=strrep(colvar,'_',' ');
COL=termtab{:,startsWith(termtab.Properties.VariableNames, colvar)};

% sort out the size and color limits/labels
d=options.sigfigs-1;
numfmt="%."+string(options.sigfigs)+"g";

szleglabs={sprintf(numfmt,min(SIZ(:))),sprintf(numfmt,max(SIZ(:)))};
if options.log_size
    switch options.log_size_mode
        case 'two'
            SIZ=log2(SIZ);
        case 'ln'
            SIZ=log(SIZ);
        case '1p'
            SIZ=log1p(SIZ);
        case 'ten'
            SIZ=log10(SIZ);
        case '-ten'
            SIZ=-log10(SIZ);
            szleglabs=fliplr(szleglabs);
        otherwise
    end
end

CBLIM=[min(COL(:)), max(COL(:))];
cbticks=linspace(CBLIM(1),CBLIM(2),options.cb_ticks);
% cbticks=[ceil(cbticks(1)*10^d)/10^d, floor(cbticks(2)*10^d)/10^d]; %make sure the ticks are in range
colleglabs=arrayfun(@(x)sprintf(numfmt,x),cbticks);
if options.log_color
    switch options.log_color_mode
        case 'two'
            COL=log2(COL);
            cbticks=log2(cbticks);
        case 'ln'
            COL=log(COL);
            cbticks=log(cbticks);
        case '1p'
            COL=log1p(COL);
            cbticks=log1p(cbticks);
        case 'ten'
            COL=log10(COL);
            cbticks=log10(cbticks);
        case '-ten'
            COL=-log10(COL);
            cbticks=-log10(cbticks);
            colleglabs=fliplr(colleglabs);
        otherwise
    end
end


% switch(sizevar)
%     case {'p_value','p','fdr'}
%         SIZ=-log10(SIZ); 
%         szleglabs=fliplr(szleglabs);
% %         sizelab="-log_{10} "
%     case {'fold_enrichment','fc','n_int'}
%         SIZ=log2(SIZ); 
% %         sizelab="log_2 " + sizelab;
% end
%
% CBLIM=[min(COL(:)),max(COL(:))];
% cbticks=[ceil(CBLIM(1)*10)/10,floor(CBLIM(2)*10)/10];
% colleglabs={sprintf(numfmt,cbticks(1)),sprintf(numfmt,cbticks(2))};
% % cbticks=[ceil(CBLIM(1)*10)/10,round(diff(CBLIM)/2*10)/10,floor(CBLIM(2)*10)/10];
% % colleglabs={sprintf(numfmt,cbticks(1)),sprintf(numfmt,cbticks(2)),sprintf(numfmt,cbticks(3))};
% switch colvar
%     case {'p_value','p','fdr'}
%         COL=-log10(COL);
%         cbticks=-log10(cbticks);
%         colleglabs=fliplr(colleglabs);
%     case {'fold_enrichment','fc','n_int'}
%         COL=log2(COL); 
%         cbticks=log2(cbticks);
% %         collab="log_2 " + collab;
% end

xlab=groupnames;
ylab=termtab.term_name;
[X,Y]=meshgrid(1:length(xlab),1:height(termtab));

if ~isempty(options.group_order)
    SIZ=SIZ(:,options.group_order);
    COL=COL(:,options.group_order);
    xlab=xlab(options.group_order);
end

% SIZ=flipud(SIZ);
% SIZ=flipud(SIZ);

% ylab - abbreviate certain terms
% ylab=strrep(ylab,'regulation of nervous','reg. of nervous');
% ylab=strrep(ylab,'regulation of neural','reg. of neural');
ylab=strrep(ylab,'Regulation','Reg.');
ylab=strrep(ylab,'regulation','reg.');
ylab=strrep(ylab,'regulating','reg.');
ylab=strrep(ylab,'positive','pos.');
ylab=strrep(ylab,'negative','neg.');
ylab=strrep(ylab,'population','pop.');
% ylab=strrep(ylab,'component','cpnt.');
ylab=strrep(ylab,'plasma membrane','PM');
ylab=strrep(ylab,'plasma-membrane','PM');
ylab=strrep(ylab,'ligand-receptor','LR');
% ylab=strrep(ylab,'Protein processing in endoplasmic reticulum','Protein processing in ER');
ylab=strrep(ylab,'Metabolism of xenobiotics by cytochrome P450','Xenobiotic metabolism - cytochrome P450');
ylab=strrep(ylab,'Extracellular matrix','ECM');
ylab=strrep(ylab,'extracellular matrix','ECM');
ylab=strrep(ylab,'Fibroblast growth factor','FGF');
ylab=strrep(ylab,'fibroblast growth factor','FGF');
ylab=strrep(ylab,'vascular endothelial growth factor','VEGF');
ylab=strrep(ylab,'platelet-derived growth factor','PDGF');
ylab=strrep(ylab,'Receptor Tyrosine Kinases','RTKs');
ylab=strrep(ylab,'Central nervous system','CNS');
ylab=strrep(ylab,'central nervous system','CNS');
ylab=strrep(ylab,'endothelial cell','EC');
ylab=strrep(ylab,'endoplasmic reticulum','ER');
% ylab=strrep(ylab,'voltage-gated calcium channel','VGCC');
% ylab=strrep(ylab,'Asparagine','Asn');

% idnum=termtab.term_id;
% idnum=strrep(idnum,'GO:','');
% idnum=strrep(idnum,'KEGG:','');
% idnum=strrep(idnum,'REAC:R-RNO-','');
% idnum=strip(idnum,'left','0');
% ylab=idnum+": "+ylab;

namelength=cellfun(@length,ylab);
longname=namelength>options.max_termlength;
abbrname=ylab(longname);
ylab(longname)=cellfun(@(x){[x(1:options.max_termlength-2),'...']},abbrname);

if ~ismember('source',termtab.Properties.VariableNames)
    options.show_source=false;
end

if options.show_source
ylab2=strrep(termtab.source,'GO:','');
% ylab2=ylab2+":"+idnum;

if ~isempty(options.ylabgroup)
    ygrp=categorical(termtab.(options.ylabgroup));
    [ugrp,ia,ic]=unique(string(ygrp));
    nygrp=length(categories(ygrp));
    if isempty(options.ylabcol)
        ycols=hsv(nygrp);
    else
        ycols=options.ylabcol;
    end
    for i=1:length(ylab2)
        ylab2(i)={sprintf('\\color[rgb]{%f,%f,%f}%s',ycols(ygrp(i),:),ylab2{i})};
        [ina,locB]=ismember(i,ia);
        if ina && options.show_ygroup_label
            ylab2(i)=cellstr(ugrp(locB)+" "+ylab2{i});
        end
    end
end
end

axpos=options.axpos;

XLIM=[1,length(xlab)]+0.5*[-1,1];
YLIM=[1,length(ylab)]+0.75*[-1,1];

if options.show_source
%dummy axis for the term source labels
ax2=axes('Position',axpos);
xlim(ax2,XLIM);
xticks(ax2,[])
yticks(ax2,1:length(ylab));
yticklabels(ax2,ylab2);
ylim(ax2,YLIM);
ax2.FontSize=8;
ax2.YAxisLocation='left';
ax2.YDir='reverse';
end

%do main axis second so it is on top
ax=axes('Position',axpos);

%%%%%%%%%%%%%bubblechart
hb=bubblechart(ax,X(:),Y(:),SIZ(:),COL(:));
bubblesize(ax,options.sizerange)
colormap(cool)

xticks(1:length(xlab))
xticklabels(xlab);
xtickangle(30);
xlim(XLIM);

yticks(1:length(ylab));
yticklabels(ylab);
ylim(YLIM);

ax.TickLabelInterpreter = 'tex';
ax.FontSize=8;
ax.YAxisLocation='right';
ax.YDir='reverse';

if ~isempty(nPerGroup)
    nPerGroup=nPerGroup(:)';
    yGridValues=cumsum(nPerGroup(1:end-1))+0.5;
    line(repmat(xlim()',1,length(yGridValues)),[1;1]*yGridValues,'color',0.5*[1,1,1],'tag','groupdiv','linewidth',0.5)
end


hcb=colorbar(ax);
hcb.Location='northoutside';
hcb.FontSize=8;
hcb.Label.String=collab; %label is a text object relative to cb axes??
hcb.Label.FontSize=8;

cbXLIM=xlim(hcb);
cbYLIM=ylim(hcb);

hcb.Location='manual';
hcb.Position(1)=axpos(1)+.0;
hcb.Position(2)=axpos(2)+axpos(4)+0.005;
hcb.Position(3)=axpos(3)-0.0;
hcb.Position(4)=0.01;
% CBLIM=hcb.Limits;
% CBLIM=sort([CBLIM,CBLIM(1)+diff(CBLIM)/2]);
% hcb.Ticks=[ceil(CBLIM(1)*10)/10,round(CBLIM(2)*10)/10,floor(CBLIM(3)*10)/10];
hcb.Ticks=cbticks;
hcb.TickLabels=colleglabs;
hcb.TickLength=0.05;

%{'';sizlab}
hbl=bubblelegend(ax,sizlab,'NumBubbles',2,'Style','horizontal');
hbl.LimitLabels=szleglabs;
hbl.Location='none';
hbl.Position(1)=axpos(1)+axpos(3)+0.075;
hbl.Position(2)=axpos(2)+axpos(4)-0.01;
hbl.Position(3)=0.05;
hbl.Position(4)=0.05;
hbl.Box='off';
hbl.FontSize=8;
hbl.Title.FontWeight='normal';

hbc.ax=ax;

if options.show_source
hbc.ax2=ax2;
end

hbc.hb=hb;
hbc.hbl=hbl;
% hbc.hblt=hblt;
hbc.cb=hcb;