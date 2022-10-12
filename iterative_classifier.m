function [type, scores, thr, M, CTs, typetop]=iterative_classifier(tcounts, genes, markers, options, markeropts)
arguments
    tcounts
    genes
    markers
    options.blockvar=ones(1,size(tcounts,2))
    options.thrsub=true(1,size(tcounts,2))
    options.thrmult=1
    options.refine=true
    options.nReps=1
    options.mintop=50
    options.fractop=.1
    options.retain_markers=false

    markeropts.summary=["min_cohen_d", "min_delta_prop"]
    markeropts.method="top";
    markeropts.ntop=10;
    markeropts.min_self_prop=0.1;
end

M=markers;
CTs=unique(M.celltype,'stable')';
nCells=size(tcounts,2);

block=categorical(options.blockvar);
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
            scores(i,blocksub)=score_genes(mnames, T, genes.name, ctrl_size=50);
        end
    end

    %get thresholds
    for i=1:length(CTs)
        this_score=scores(i,:);
        thr(i,1)=geneThresholdOtsu(this_score(options.thrsub))*options.thrmult;
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

        ges=GroupedExpressionSummary(genes,tcounts,typetop,block);
        ges.computeEffectSizes();


        M=table;
        for i=1:length(markeropts.summary)
            M=[M;ges.topMarkers(markeropts.summary(i),markeropts.method, CTs, ...
                nTop=markeropts.ntop, min_self_prop=markeropts.min_self_prop)];
%             M=[M;ges.topMarkers(markeropts.summary(i),markeropts.method,, ...
%                 nTop=markeropts.ntop, min_self_prop=markeropts.min_self_prop)];
        end
        M=M(:,1:2);
        if options.retain_markers
            M=[M;Mlast];
%             M=[M;markers(ismember(markers.celltype,CTs),:)];
        end
        M=unique(M,'rows','stable');
        M=sortrows(M,["celltype","gene"]);

%         equalM=isequal(Mlast,M)
        summary(categorical(M.celltype,CTs,CTs)')

        Mlast=M;
    end
end