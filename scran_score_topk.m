function result = scran_score_topk(scorestruct, nTop, selectby, options)
arguments
    scorestruct
    nTop
    selectby="min.logFC.detected"
    options.selftypes = {}
    options.combineby="union"
    options.self_minprop=0
    options.other_maxprop=1
    options.direction="up"
end

selftypes=options.selftypes;
if isempty(selftypes)
    selftypes=scorestruct.pairings;
    options.selftypes=selftypes;
end

minlFC=0;
minAUC=0.5;

result=options;
result.k=nTop;
result.selectby=selectby;

for i=1:length(selftypes)
    thisCT=scorestruct.(selftypes{i});
    %prefilter
    filt=thisCT.("self.detected")>=options.self_minprop & thisCT.("other.detected")<=options.other_maxprop;
    thisCT=thisCT(filt,:);
    clear thisres stats
    IX={};
    for j=1:length(selectby)
        [m,ix]=maxk(thisCT.(selectby(j)),nTop);
        if contains(selectby(j),"logFC") && any(m<0)
            warning('negative logFC')
            ix(m<0)=[];
        end
        if contains(selectby(j),"AUC") && any(m<0.5)
            warning('AUC<0.5')
            ix(m<0)=[];
        end
        IX{j}=ix;
        thislab=strrep(selectby(j),".","_");
        thisres.(thislab)=thisCT(ix,:);
        stats.minES(j,1)=min(m);
        stats.maxES(j,1)=max(m);
        stats.minSelfProp(j,1)=min(thisCT{ix,"self.detected"});
        stats.maxSelfProp(j,1)=max(thisCT{ix,"self.detected"});
        stats.minOtherProp(j,1)=min(thisCT{ix,"other.detected"});
        stats.maxOtherProp(j,1)=max(thisCT{ix,"other.detected"});
    end
    uix=unique(cat(1,IX{:}));
    switch options.combineby
        case "union"
            ix=uix;
        case "intersect"
            ix=uix;
            for j=1:size(IX,2)
                ix=intersect(ix,IX{j});
            end
    end
    selected=thisCT(ix,:);

    thisres.stats=struct2table(stats,"RowNames",selectby);
    thisres.summary=sortrows(selected,'gene'); %alphabetical
    thisres.genes=thisres.summary.gene;
    result.(selftypes{i})=thisres;
end