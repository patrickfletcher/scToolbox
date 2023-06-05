function result=iterative_Otsu_classifier(tcounts, genes, markers, options)
arguments
    tcounts
    genes
    markers
    options.blockvar=[]
    options.ctrl_genes=50
    options.thrmult=1
    options.refine=true
    options.nReps=1
    options.mintop=50
    options.fractop=.1
    options.retain_markers=false

    options.ref_subset=[]

    options.summary="min_cohen_d"
    options.summary_min=0.5
    options.method="top"
    options.ntop=10;
    options.thr=5;
    options.p=95;
    options.min_self_prop=0
    options.max_other_prop=1
end
% Use gene scores + thresholding to call each cell type - a cell can be
% assigned 0, 1, or many cell types.

%what about using score ranks? generate ecdf: then scores are all in [0,1]

M=markers;
CTs=unique(M.celltype,'stable')';
CTo=[CTs,"Unc"];
% CTo=CTs;
nCells=size(tcounts,2);

block=options.blockvar;
if isempty(block)
    block=true(size(tcounts,2),1);
end
block=categorical(block);
blocks=categories(block);


ref_subset=options.ref_subset;
if isempty(ref_subset)
    ref_subset=true(size(tcounts,2),1);
end
refix=find(ref_subset);

if options.refine
    options.nReps=max(2,options.nReps);
end

typelast=strings(nCells,1);
Mlast=M(:,["celltype","gene"]);
M0=M(:,["celltype","gene"]);

% for i=1:length(options.summary)
%     Mlast.(options.summary(i))=zeros(size(M.gene));
% end

for r=1:options.nReps

    % get marker scores
    scores=zeros(length(CTs),nCells);
    scorethr=false(length(CTs),nCells);
    thr=zeros(length(CTs),length(blocks));
    for j=1:length(blocks)
        blocksub=block==blocks{j};
        T=tcounts(:,blocksub);
        for i=1:length(CTs)
            thisCT=CTs(i);
            mnames=M.gene(M.celltype==thisCT);
            this_score=score_genes(mnames, T, genes.name, ctrl_size=options.ctrl_genes);
            scores(i,blocksub)=this_score;

            % Identify top scoring cells per marker set
            % - Otsu thresholds
            thr(i,j)=geneThresholdOtsu(this_score)*options.thrmult;
            scorethr(i,blocksub)=this_score>thr(i,1);
        end
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

    if isequal(typelast,type)
        disp("finished in " + string(r) + " iterations.")
        break
    end
    typelast=type;

    % refine markers
    if options.refine && r<options.nReps
        
        TOP=[];
        %fraction applied to each type, or total cells?
%         Q = quantile(scores,1-options.fractop,2)
        for i=1:length(CTs)
            thistype=type==CTs{i} & ref_subset;
            keepn=max(options.mintop, round(options.fractop*nnz(thistype)));
            keepn=min(keepn, nnz(thistype));
            [~, ix]=maxk(scores(i,ref_subset),keepn); %may include amb...
            topix=refix(ix);
            % [~, ix]=maxk(scores(i,thistype),keepn); %exclude amb...
            % subix=find(thistype);
            % topix=subix(ix);
            TOP=[TOP;topix(:)];
        end
        TOP=unique(TOP);
        discard=setdiff(1:length(type),TOP);
        
        typetop=removecats(type,["Amb","Unc"]);
        typetop(discard)=missing();
        unc=find(type=="Unc"); 
        keepn=max(options.mintop, round(options.fractop*nnz(type=="Unc")));
        kix=randi(length(unc),keepn,1);
        typetop(unc(kix))="Unc";
        % typetop(type=="Unc")="Unc";
        typetop=removecats(typetop);
        % summary(typetop)

        % get the new marker list based on top scoring cells
        ges=GroupedExpressionSummary(genes,tcounts,typetop,block=block);
        ges.computeEffectSizes();

        M=table;
        for i=1:length(options.summary)
            M=[M;ges.topMarkers(options.summary(i),options.method, CTs, CTo, ...
                  nTop=options.ntop, summary_min=options.summary_min(i),...
                  min_self_prop=options.min_self_prop, max_other_prop=options.max_other_prop)];
        end
        M=M(:,["celltype","gene"]);
        % M=M(:,["celltype","gene",options.summary]);
        if options.retain_markers
            M=[M;M0];
            % M=[M;Mlast];
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