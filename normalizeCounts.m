function [ncounts,sfs]=normalizeCounts(counts, options)
arguments
    counts
    options.scale = []
    options.center_measure = "median"
    options.sizefactors=[]
end
%normalize counts

%TODO
% - batch handling: multibatchnorm?
% - exclude genes making up >max_frac of a cell's expression? scanpy [e.g. poor quality cells/RBCs?]
% - subset of genes used to compute trend?

% could support others... https://en.wikipedia.org/wiki/Central_tendency 
switch options.center_measure
    case "mean"
        center_measure = @mean;
    case "median"
        center_measure = @median;
end

%sum counts per cell
counts_per_cell=full(sum(counts,1));

% rescale the normalized counts to a) original scale (measure?) or b) a specified target scale
scale=options.scale;
if isempty(scale) %sf_center_target = 1 gives original scale using center_measure
    scale = center_measure(counts_per_cell);
end

%basic size factor: divide by counts_per_cell --> fraction of cell total counts per gene
sfs=options.sizefactors;
if isempty(sfs)
    sfs=counts_per_cell/scale;
else
    % adjust user sfs to get correct scale when using "center_measure"
    % - for original data scale: sf_center_target = 1 (center the sfs at 1)
    % - for specified scale: sf_center_target = center_measure(counts_per_cell)/scale
    sf_center_target = center_measure(counts_per_cell)/scale;
    sf_center = center_measure(sfs);
    sfs = sfs/sf_center * sf_center_target;
end

sfs=sfs(:)';
if issparse(counts)
%     A=spfun(@(x)1./x,sfs);
    [i,j,a]=find(counts);
    [m,n]=size(counts);
    an=zeros(size(a));
    for k=1:length(a)
        an(k)=a(k)./sfs(j(k));
    end
    ncounts=sparse(i,j,an,m,n);
else
    ncounts=counts./sfs;
end

% verify
% scalef = center_measure(full(sum(ncounts,1)));
% [scale,scalef]
% (scalef-scale)/scale*100 %relative error




%???

    % options.min_mean = 0.1
    % options.batch = []
    % options.max_frac=[]

% max_frac=options.max_frac;
% if ~isempty(max_frac)
%     %redo the counts_per_cell
%     high_frac=full(sum(counts>max_frac*counts_per_cell,2));
%     counts_per_cell=full(sum(counts(~high_frac,:),1));
%     scale=median(counts_per_cell);
% end

% batch = options.batch;
% if isempty(batch)
%     batch = ones(1,size(counts,2));
% end
% batch = findgroups(batch);

% gsub = true(size(counts,1),1);
% if ~isempty(options.min_mean)
%     gsub = 
% end

% scaleb = splitapply(@median, counts_per_cell(:), batch(:));

    % %should do this per batch? find the batch with mean sfs=1?
    % sfs_b = splitapply(@mean, sfs(:), batch(:));
    % [~,ix]=min(abs(sfs_b-1));
    % old_scale = scale_b(ix);
    % sfs = sfs.*old_scale./scale;
    
% %batch: multibatchnorm.
% % 1. compute sfs=median(counts_per_cell) for all. (libSizefactors - above)
% % 2. for each batch: sf / mean(sf), mean(counts/sf)
% % 3. using only genes mean_ncounts>min_mean, find smallest sf. Shrink sizefactors of all batches to match the smallest.
% if ~isempty(options.batch)
%     b = findgroups(options.batch);
%     sfs_bm = splitapply(@mean, sfs, batch);
%     sfs_ratio = splitapply(@(x){x./mean(x)},sfs,batch); 
%     sfs_ratio=cat(2,sfs_ratio{:});
% end