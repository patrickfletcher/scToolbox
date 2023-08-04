function hvqc=plot_violin_QC(cells, qc_data, blockvar, opts)
arguments
    cells
    qc_data
    blockvar = []
    opts.varlabs = [] 
    opts.gcols = []
    opts.fig = []
end

varnames=[qc_data(:).name]; 
do_log=[qc_data(:).do_log];

if isempty(blockvar)
    blockvar=ones(size(cells,1), 1);
end
varlabs=opts.varlabs;
if isempty(varlabs)
    varlabs=varnames;
end

data=cells{:,varnames};
data(:,do_log)=log10(data(:,do_log));

gcols=opts.gcols;
if isempty(gcols)
    gcols = turbo(length(unique(blockvar)));
end

hvm=violinmatrix(data, blockvar, varlabs, gcols, ...
    grid='on',ShowMedian='on', ShowBox='on', ...
    ShowData='off', fig=opts.fig);

for i=1:length(hvm.v(:))
    hvm.v(i).ViolinAlpha=0.9;
    hvm.v(i).MedianPlot.SizeData=5;
    hvm.v(i).ScatterPlot.Marker='x';
end

K=length(unique(blockvar));
for i=1:length(varnames)
    x=[0.66;1.33]+(1:K)-1;
    yup=repmat(qc_data(i).stats.upper(:)',2,1);
    ylo=repmat(qc_data(i).stats.lower(:)',2,1);
    if qc_data(i).do_log
%         ax(i).YScale='log';
        yup=log10(yup);
        ylo=log10(ylo);
    end
    if qc_data(i).type=="upper" || qc_data(i).type=="both" 
        hlhi{i}=line(hvm.ax(i),x,yup,'color',0.6*[1,1,1]);
    end
    if qc_data(i).type=="lower" || qc_data(i).type=="both" 
        hllo{i}=line(hvm.ax(i),x,ylo,'color',0.6*[1,1,1]);
    end
%     line(ax(i),[min(x(:)),max(x(:))],median(data(:,i))*[1,1],'color',0.7*[1,1,1])
    ylim(hvm.ax(i),'padded')
end

hvqc = hvm;
hvqc.hlhi=hlhi;
hvqc.hllo=hllo;
