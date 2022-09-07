function result = classify_cells(genes, tcounts, initids, initmarkers, options)
arguments
    genes
    tcounts
    initids
    initmarkers = []
    options.marker_selectby = "min_cohen_d"
    options.marker_method = "top"
    options.marker_nTop = 20
    options.refine = true
    options.refine_ptile = 95
end

CTs=categories(initids);
nCells=size(tcounts,2);

if isempty(initmarkers)
    ges=GroupedExpressionSummary(genes,tcounts,initids);
    ges.computeEffectSizes();
    
    M=table;
    for i=1:length(options.marker_selectby)
        M=ges.topMarkers(options.marker_selectby(i), options.marker_method, nTop=options.marker_nTop);
    end
    M=M(:,1:2);
    M=unique(M,'rows');
end

% compute ct-scores 
scores=zeros(length(CTs),nCells);
for i=1:length(CTs)
    thisCT=CTs(i);
    mnames=M.gene(M.celltype==thisCT);
    scores(i,:)=score_genes(mnames,tcounts, genes.name, ctrl_size=50);
end

% select thresholds from score distributions, assign labels
% - otsu, quantile, ...
% - SingleR uses correlation values as scores, allowing numeric
%   first/second comparisons. 
% - Problem with genescores is they are not on the same scale (or are
%   they?) - normalize them???
scorethr=false(length(CTs),nCells);
thr=zeros(length(CTs),1);
for i=1:length(CTs)
    this_score=scores(i,:);
    thr(i,1)=geneThresholdOtsu(this_score);
    scorethr(i,:)=this_score>thr(i,1);
end

utype=sum(scorethr,1)==1;
type=strings(nCells,1);
CTrep=repmat(CTs,1,nCells);
uCTs=CTrep(scorethr(:,utype));
type(utype)=uCTs;
type(sum(scorethr,1)==0)="Unc";
type(sum(scorethr,1)>1)="Amb";

ctplus=[CTs,"Unc","Amb"];
type=categorical(type,ctplus,ctplus);

summary(type')

if options.refine
    ges=GroupedExpressionSummary(genes,tcounts,type);
    ges.computeEffectSizes();
    
    M=table;
    for i=1:length(options.marker_selectby)
        M=ges.topMarkers(options.marker_selectby(i), options.marker_method, nTop=options.marker_nTop);
    end
    M=M(:,1:2);
    M=unique(M,'rows');
end
