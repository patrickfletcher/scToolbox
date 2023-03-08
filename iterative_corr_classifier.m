function result=iterative_corr_classifier(tcounts, genes, markers, options)
arguments
    tcounts
    genes
    markers
    options.blockvar=[]

    options.mintop=50
    options.fractop=.1
    options.ctrl_genes=50

    options.quantile=0.8
    options.min_delta = 0.1

    options.refine=true
    options.nReps=1
    options.retain_markers=false

    options.summary=["min_cohen_d", "min_delta_prop"]
    options.method="top";
    options.ntop=10;
    options.min_self_prop=0.1;

    options.debug_output=false;
end
% SingleR approach
% - top scoring cells to define "reference" cells per type
% -- pooling, vector quantized, sqrt(n) clusters (kmeans) to get multiple ref samples per type
% - compute the Spearman correlation between each cell's expression profile and
% each cell type reference. This gives a correlation coefficient to
% each cell type for each cell (also p-val). scores based on this.
%
% scores:
% - per-label score is a fixed quantile (by default, 0.8) of the correlations across all samples with that label
%
% diagnostics:
% - delta - difference between a cell's score for the assigned label and the median score across all labels
% - 
%
% pruning: minimum delta? minimum next delta?

% once types assigned, can switch from topN to fracTop for reference cells?

M=markers;
CTs=unique(M.celltype,'stable')';
nCells=size(tcounts,2);

block=options.blockvar;
if isempty(block)
    block=ones(1,size(tcounts,2));
end
block=categorical(block);
blocks=categories(block);

if options.refine
    options.nReps=max(2,options.nReps);
end

typelast=strings(nCells,1);
Mlast=M;

type = typelast;

for r=1:options.nReps
    allMarkers = M.gene;
    allMarkerIx = getGeneIndices(allMarkers,genes.name);

    % get marker gene scores: reference cell set via marker gene scores
    mscores=zeros(length(CTs),nCells);
    for j=1:length(blocks)
        blocksub=block==blocks{j};
        T=tcounts(:,blocksub);
        for i=1:length(CTs)
            thisCT=CTs(i);
            mnames=M.gene(M.celltype==thisCT);
            mscores(i,blocksub)=score_genes(mnames, T, genes.name, ctrl_size=options.ctrl_genes);
        end
    end
    
    % Identify top scoring cells per marker set: keep a fixed N?
    % - need to verify that these don't include Amb
    RefIx={};
    alltopix=[];
    for i=1:length(CTs)
        %first pass: type is not populated, so will select n=mintop
        keepn=max(options.mintop, round(options.fractop*nnz(type==CTs{i})));
        [~, topix]=maxk(mscores(i,:),keepn,2);
        multiref=ismember(topix,alltopix);
%         any(multiref)
        topix = topix(~multiref);
        RefIx{i}=topix;
        alltopix=[alltopix,topix];
    end
    
    if options.debug_output
        result.RefIx=RefIx;
        result.alltopix=alltopix;
    end

    % get the corrs (all cells corr to ref cells)
    % - use full set of marker genes for correlation
    scores=zeros(length(CTs),nCells);
    T=tcounts(allMarkerIx,:); %cells as rows
    allZeros=all(T==0,1);
    for i=1:length(CTs)
        Tref=tcounts(allMarkerIx,RefIx{i}); %cells as rows
%         allZref=all(Tref==0,1); nnz(allZref)
        [RHO, P]=corr(T, Tref, type="Pearson", tail='right');
%         [RHO, P]=corr(T, Tref, type="Spearman", tail='right');
        % rho of Ref with Ref should be high, but may not be...
%         RHO(allZeros,:)=0; %cells not expressing the markers get all zeros
%         scores(i,:)=prctile(RHO, 100, 2);
        scores(i,:)=quantile(RHO, options.quantile,2); %how many ref cells do we require high correlation with?
        if options.debug_output
            result.RHO{i}=RHO;
            result.P{i}=P;
        end
    end

    [maxscores,mix]=maxk(scores,1,1);

    type=categorical(CTs(mix),CTs,CTs);
    summary(type')
    typeRef = type;
    typeRef(setdiff(1:nCells, alltopix))=missing();
    summary(typeRef')

    medscores=median(scores,1);
    deltas=scores-medscores;

    [maxdeltas,mdix]=maxk(deltas,2,1);
    diffMaxDeltas = diff(flipud(maxdeltas),1,1);

    % SingleR: small delta outliers per type = unc - this relies on distribution across cells
    % ambiguous: not Unc and top two deltas are too close
%     [~, qcdat]=is_outlier(deltas', "lower", type, nmads=2, clip_min=-Inf);

%     thisAmb = diffMaxDeltas < options.min_delta;
%     type(thisAmb) = "Amb";

    for i=1:length(CTs)
        thisType = type==CTs{i};

%         thisAmb = thisType & diffMaxDeltas < options.min_delta;
%         type(thisAmb) = "Amb";

%         outlier=isoutlier(deltas(i,thisType),"median");
        outlier=is_outlier(deltas(i,thisType)', "lower", nmads=2, clip_min=-Inf);
        thisUnc = find(thisType(outlier));
%         thisUnc = thisType & qcdat(i).outlier;
%         thisUnc = thisType & maxdeltas(1,:) < options.min_delta;
        type(thisUnc) = "Unc";
    end

    type(allZeros)="Unc";
    
    summary(type')

    % resolve Unc/Amb (SingleR pruning)
    % - threshold on maxdelta? maxdelta < 0.1 = Unc?
    % - more than one maxdelta above threshold: Amb?

    deltathr = deltas > options.min_delta;
    deltathr2 = maxdeltas(1,:) > maxdeltas(2,:) + options.min_delta;
    unccond = deltathr & deltathr2;
%     ambcond = deltathr & 
    utype=sum(deltathr & deltathr2,1)==1;
    type=strings(nCells,1);
    CTrep=repmat(CTs,1,nCells);
    uCTs=CTrep(deltathr(:,utype));
    type(utype)=uCTs;
    type(sum(deltathr,1)==0)="Unc";
    type(sum(deltathr,1)>1)="Amb";
    
    ctplus=[CTs,"Unc","Amb"];
    type=categorical(type,ctplus,ctplus);

    summary(type')

    if isequal(typelast,type)
        break
    end
    typelast=type;
% 
%     % refine markers
%     if options.refine && r<options.nReps
%         
%         typetop=removecats(type,["Amb","Unc"]);
%         typetop(~refcells)=missing(); 
%         unc=find(type=="Unc"); 
%         keepn=max(options.mintop, round(options.fractop*nnz(type=="Unc")));
%         kix=randi(length(unc),keepn,1);
%         typetop(unc(kix))="Unc";
%         typetop=removecats(typetop);
% 
%         % get the new marker list based on top scoring cells
%         ges=GroupedExpressionSummary(genes,tcounts,typetop,block);
%         ges.computeEffectSizes();
% 
%         M=table;
%         for i=1:length(options.summary)
%             M=[M;ges.topMarkers(options.summary(i),options.method, CTs, CTs, ...
%                   nTop=options.ntop, min_self_prop=options.min_self_prop)];
%         end
%         M=M(:,1:2);
%         if options.retain_markers
%             M=[M;Mlast];
% %             M=[M;markers(ismember(markers.celltype,CTs),:)];
%         end
%         M=unique(M,'rows','stable');
%         M=sortrows(M,["celltype","gene"]);
% 
% %         equalM=isequal(Mlast,M)
%         summary(categorical(M.celltype,CTs,CTs)')
% 
%         Mlast=M;
%     end
end

result = options;
result.type = type;
result.marker_scores = mscores;
result.scores = scores;
result.medscores = medscores;
result.deltas = deltas;
result.markers = M;
result.CTs = CTs;