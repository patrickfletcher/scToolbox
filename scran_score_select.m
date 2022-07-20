function result = scran_score_select(scorestruct, options)
arguments
    scorestruct
    options.selfnames={}
    options.othernames={}
    options.self_minprop=0.25
    options.other_maxprop=1
    options.effect_type="min.AUC"
    options.select_method="value"
    options.threshold=1
    options.direction="up"
    options.pairwise_logic="any"
    options.sortby="rank.logFC.detected"
    options.sortdir="ascend"
end
% scorestruct is result of "scran_score_markers"
% NOTE: the summaries from scoreMarkers don't true combining of effects at
% pairwise level:
% - min.AUC>thr & min.logFC>thr NOT= all(AUC>thr & logFC>thr, 2)

% select_methods:
% "value" - select using threshold on ES values
% "topk" - select larges K values of ES  *** this is simple externally

% - maybe: if no thresh, set defaults for up/down (AUC=0.5, logFC=0)
% - filter on thresholds, then take topk of remaining if specified?

%TODO: more convenient to just return simpler results? selected ES, mean,
%prop. can always add to scorestruct externally...

%TODO: self/other prop threshold should use min/max on prop table
% - scoreMarkers' other.detected column is mean of others...
% - even the prop/mean tables are corrected, so may differ from raw version
% of same calculation in individual blocks.

%TODO: support fullstats usage - custom pairings/combinations 
% - when fullstats+self/other, collect full ES, apply combine rule...
% - specify combo of effect_type + summary_type (min/median/.../full??)
nCrit=length(options.effect_type);

% result=scorestruct;

%TODO: handle specifying self/other names to do specific pairings.
selftypes=options.selfnames;
if isempty(selftypes)
    selftypes=scorestruct.pairings;
end

if nCrit>1 && length(options.threshold)==1
    options.threshold=repmat(options.threshold,nCrit,1);
end
if nCrit>1 && length(options.direction)==1
    options.direction=repmat(options.direction,nCrit,1);
end

result.selectopts=options;
for i=1:length(selftypes)
    thisCT=scorestruct.(selftypes{i});
    thisCT=sortrows(thisCT,options.sortby,options.sortdir);
    otherprop=scorestruct.prop{:,setdiff(selftypes,selftypes{i})};
    
    pw_filt=false(size(thisCT.gene,1),nCrit);
    for j=1:nCrit
        pw_up=thisCT.(options.effect_type(j))>options.threshold(j);
        pw_down=thisCT.(options.effect_type(j))<options.threshold(j);
        switch options.direction{j}
            case "up"
                pw_filt(:,j)=pw_up(:);
            case "down"
                pw_filt(:,j)=pw_down(:);
            case "both"
                pw_filt(:,j)=pw_up(:)|pw_down(:);
        end
    end

    switch options.pairwise_logic
        case "any"
            comb_pw_filt=any(pw_filt,2);
        case "all"
            comb_pw_filt=all(pw_filt,2);
    end

    self_filt=thisCT.("self.detected")>options.self_minprop;
    other_filt=all(otherprop<options.other_maxprop,2); %ideally would be simultaneously pairwise with fullstats...

    filt = self_filt(:) & other_filt(:) & comb_pw_filt(:);

    filttab=table;
    filttab.("self.detected")=self_filt(:);
    filttab.("other.detected")=other_filt(:);
    filttab=[filttab,array2table(pw_filt,VariableNames=options.effect_type)];
    filttab.("combined.pw")=comb_pw_filt(:);
    filttab.Properties.RowNames=thisCT.gene;
    filttab.finalfilt=filt;

    result.(selftypes{i})=thisCT(filt,:);
    result.(selftypes{i}+"_filters")=filttab;
end