function result=iterative_cdf_classifier(tcounts, genes, markers, options)
arguments
    tcounts
    genes
    markers
    options.blockvar=[]
    options.ctrl_genes=50
    options.scorethr=0.8
    options.refine=true
    options.nReps=1
    options.mintop=50
    options.retain_markers=true

    options.ref_subset=[]

    options.summary=["min_cohen_d", "min_delta_prop"]
    options.method="top";
    options.ntop=10;
    options.min_self_prop=0.1;
end
% Marker gene scores converted to ecdf: then scores are all in [0,1]

M=markers;
CTs=unique(M.celltype,'stable')';
CTo=[CTs;"Unc"];
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
subix=find(ref_subset);

if options.refine
    options.nReps=max(2,options.nReps);
end

typelast=strings(nCells,1);
Mlast=M;

for r=1:options.nReps

    % get marker scores
    scores=zeros(length(CTs),nCells);
    scorethr=false(length(CTs),nCells);
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
    % - Otsu thresholds
    for i=1:length(CTs)
        this_score=ecdf(scores(i,:));
        scorethr(i,:)=this_score>options.scorethr;
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
        %fraction applied to each type, or total cells?
        for i=1:length(CTs)
            thistype=type==CTs{i} & ref_subset;
            keepn=max(options.mintop, round(options.fractop*nnz(thistype)));
            [~, topix]=maxk(scores(i,ref_subset),keepn,2);
            TOP=[TOP,topix];
        end
        TOP=unique(TOP);
        keep=ismember(1:length(type),TOP);
        
        %add Unc as a group to compare CTs to?
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
            M=[M;ges.topMarkers(options.summary(i),options.method, CTs, CTo, ...
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