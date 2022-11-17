function hhm = heatmap_plot(cvals, opts)
arguments
    cvals
    opts.xorder=[]
    opts.yorder=[]
    
    opts.cmap=[]
    opts.isdiverging=false

    opts.title=[] %label to use as title
    opts.fig=[] %specify figure to plot in (if empty, new figure)
    opts.margins=[0.1,0.1,0.1,0.1] %[left, bottom, right, top]

    opts.cbgap=0.01 %space between ax and cb
    opts.cbdims=[0.3,0.01] %[length, width]
    opts.cblabel=[]
    opts.cbLoc {mustBeMember(opts.cbLoc,["east","west","north","south"])}='east'
    opts.cbJust {mustBeMember(opts.cbJust,["low","midlo","mid","midhi","high"])}='mid'
    opts.cbDigits=1
end

%indicate ids using scatter with  'Clipping', 'off'?

l=opts.margins(1);
b=opts.margins(1);
r=opts.margins(1);
t=opts.margins(1);

ax=tight_subplot(1,1,1,0,[b,t],[l,r]);

% need an attached mini-axis for cell annotation
% ax=tight_subplot(1,1,1,0,[b,t],[l,r]);
% imagesc(cg.clust.clusterID(1:10:end,1)')
% colormap(opts.gcols)

IM=vals;
IM=smoothdata(IM,2,"loess",100);
IM(IM<0)=0;

IM=normalize(IM,2,'range');

% clustering genes
IMs=smoothdata(IM,2,"loess",500);
D=pdist(IMs,'euclidean');
Z=linkage(D,"ward");
o=optimalleaforder(Z,D,Criteria="group");

%order by peak
% [~,pix]=max(IMs,[],2);
% [~,o]=sort(pts(pix));

IM=IM(o,:);
glabs="\it"+glist(o);


imagesc(ax,pts,1:length(glist),IM)

% axis("tight")
yticks(1:length(glabs))
yticklabels(glabs)

if isdiverging
    cmap=split_cmap;
else
    cmap=cbrewer('seq','Reds',256); cmap(1,:)=[1,1,1];%cmap(1,:)=rgb('lightgray');
end
colormap(cmap)