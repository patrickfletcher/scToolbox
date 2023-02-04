function result=iterative_corr_classifier(tcounts, genes, markers, options)
arguments
    tcounts
    genes
    markers
    options.blockvar=[]
    options.ctrl_genes=50
    options.refine=true
    options.nReps=1
    options.mintop=50
    options.fractop=.1
    options.retain_markers=false

    options.summary=["min_cohen_d", "min_delta_prop"]
    options.method="top";
    options.ntop=10;
    options.min_self_prop=0.1;
end
% SingleR approach
% - pool top scoring cells to define a "reference" cell
% - compute the Spearman correlation between its expression profile and
% that of each cell type reference This gives a correlation coefficient to
% each cell type for each cell (also p-val). Assign to 0, 1, or more cell
% types based on that.
% - need a way to initialize references: marker gene score topk?
%
% scores:
% - per-label score is a fixed quantile (by default, 0.8) of the correlations across all samples with that label
% - benefit: scores are comparable across cell types (in [0,1]).
%
% diagnostics:
% - deltas - difference between the score for the assigned label and the median across all labels for each cell


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

for r=1:options.nReps

    % get marker scores
    scores=zeros(length(CTs),nCells);
    scorethr=false(length(CTs),nCells);
    thr=zeros(length(CTs),1);
    for j=1:length(blocks)
        blocksub=block==blocks{j};
        T=tcounts(:,blocksub);
        for i=1:length(CTs)
            thisCT=CTs(i);
            mnames=M.gene(M.celltype==thisCT);
            scores(i,blocksub)=score_genes(mnames, T, genes.name, ctrl_size=options.ctrl_genes);
        end
    end

    % Identify top scoring cells per marker set

    
    for i=1:length(CTs)
        this_score=scores(i,:);
        thr(i,1)=geneThresholdOtsu(this_score)*options.thrmult;
        scorethr(i,:)=this_score>thr(i,1);
    end

    % - alt: using correlation to reference built from topK scoring cells
    
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

    if isequal(typelast,type)
        break
    end
    typelast=type;

    % refine markers
    if options.refine && r<options.nReps
        
        TOP=[];
        for i=1:length(CTs)
            keepn=max(options.mintop, round(options.fractop*nnz(type==CTs{i})));
            [~, topix]=maxk(scores(i,:),keepn,2);
            TOP=[TOP,topix];
        end
        TOP=unique(TOP);
        keep=ismember(1:length(type),TOP);
        
        typetop=removecats(type,["Amb","Unc"]);
        typetop(~keep)=missing(); 
        unc=find(type=="Unc"); 
        keepn=max(options.mintop, round(options.fractop*nnz(type=="Unc")));
        kix=randi(length(unc),keepn,1);
        typetop(unc(kix))="Unc";
        typetop=removecats(typetop);

        % get the new marker list based on top scoring cells
        ges=GroupedExpressionSummary(genes,tcounts,typetop,block);
        ges.computeEffectSizes();

        M=table;
        for i=1:length(options.summary)
            M=[M;ges.topMarkers(options.summary(i),options.method, CTs, CTs, ...
                  nTop=options.ntop, min_self_prop=options.min_self_prop)];
        end
        M=M(:,1:2);
        if options.retain_markers
            M=[M;Mlast];
%             M=[M;markers(ismember(markers.celltype,CTs),:)];
        end
        M=unique(M,'rows','stable');
        M=sortrows(M,["celltype","gene"]);

%         multimarkers? genes listed for more than one cell type...
%         multimarkers=length(M.gene)-length(unique(M.gene));

%         equalM=isequal(Mlast,M)
        summary(categorical(M.celltype,CTs,CTs)')

        Mlast=M;
    end
end

result = options;
result.type = type;
result.scores = scores;
result.thr = thr;
result.markers = M;
result.CTs = CTs;
result.typetop = typetop;
result.refs = refs;