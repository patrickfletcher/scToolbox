function [outliers, result] = find_outliers(qcvals, tf, ctnames, params)
% arguments
% end

% general function to find outliers for cell-wise metric (qcvals), per cell
% type. Cell types are defined by tf - true/false table with
% rows=condition, cols=cells. ctnames required for manual adjustment of
% thresholds.

%use cases: 
% - tf+ctnames - matrix of celltype conditions, ctname per row - allows
% ambiguous
% - cellID categorical (no ctnames needed, extract). mutually exclusive
% groups

result.params=params; %store params and results together idiom

% individual type outlier cutoffs
ctcounts=sum(tf,2);
outsub=false(size(tf));
outliers=false(size(qcvals));
for i=1:size(tf,1)
    thissub=tf(i,:);
    thisctname=ctnames{i};
    thisqc=qcvals(thissub);
    if params.logval
        thisqc=log10(thisqc);
    end
    
    %stats (for export too)
    means(i)=mean(thisqc);
    medians(i)=median(thisqc);
    stds(i)=std(thisqc);
    mad0s(i)=mad(thisqc,0);
    mad1s(i)=mad(thisqc,0);
    
    lozwarning = warning('off', 'MATLAB:log:logOfZero');
    geomeans(i)=exp(mean(log(thisqc))); %equiv to geomean
    geostds(i)=exp(std(log(thisqc)));
    warning(lozwarning);
            
    %make manual settings apply to the threshold parameter, not fixed val?
    params=result.params;
    if ismember(thisctname,params.manual_low.names)
        params.ndev_lo=params.manual_low.vals(params.manual_low.names==thisctname);
    end
    if ismember(thisctname,params.manual_high.names)
        params.ndev_hi=params.manual_high.vals(params.manual_high.names==thisctname);
    end
    
    switch params.method
        case 'fixed'
            lowthr(i)=params.lowthr;
            hithr(i)=params.hithr;
            
        case 'meanstd'
            centerqc(i)=mean(thisqc);
            spreadqc(i)=std(thisqc);
            lowthr(i)=centerqc(i)-params.ndev_lo*spreadqc(i); 
            hithr(i)=centerqc(i)+params.ndev_hi*spreadqc(i);  
            
        case 'geomeanstd'
            lozwarning = warning('off', 'MATLAB:log:logOfZero');
            centerqc(i)=mean(log(thisqc)); %equiv to geomean
            spreadqc(i)=std(log(thisqc));
            warning(lozwarning);
            lowthr(i)=exp(centerqc(i)-params.ndev_lo*spreadqc(i)); 
            hithr(i)=exp(centerqc(i)+params.ndev_hi*spreadqc(i)); 
            
        case 'medmad'
            centerqc(i)=median(thisqc);
            spreadqc(i)=mad(thisqc,params.madflag);
            lowthr(i)=centerqc(i)-params.ndev_lo*spreadqc(i); 
            hithr(i)=centerqc(i)+params.ndev_hi*spreadqc(i); 
            
        case 'prctile'
            lowthr(i)=prctile(thisqc,params.lowptile);
            hithr(i)=prctile(thisqc,100-params.hiptile);
            
    end
    
    %don't impose thresholds (except clipvals) on low-pop cts
    if isempty(thisqc) || length(thisqc)<params.min_cells
        lowthr(i)=params.lowthr_clipval;
        hithr(i)=params.hithr_clipval;
        continue
    end
    
    lowthr(i)=max(lowthr(i),min(thisqc)); %clamp thresholds - don't go beyond data max
    hithr(i)=min(hithr(i),max(thisqc)); 
    
    lowthr(i)=max(lowthr(i),params.lowthr_clipval); %clamp to manual clip-vals thresholds
    hithr(i)=min(hithr(i),params.hithr_clipval); 
    
%     if ismember(thisctname,params.manual_low.names)
%         lowthr(i)=params.manual_low.vals(params.manual_low.names==thisctname);
%     end
%     if ismember(thisctname,params.manual_high.names)
%         hithr(i)=params.manual_high.vals(params.manual_high.names==thisctname);
%     end
    
    %special value for "Unc": lowthr=min(others), hithr=max(others)
    % Assumes Unc is last. 
    % alt versions? 
    % - cell count weigthed thr
    % - min/max of kept values observed in any other type
    if ctnames(i)=="Unc"
        if isempty(params.unctypes)
            otherix=ctnames~="Unc"; %use all others for weighted thresh
        else
            otherix=ismember(ctnames,params.unctypes);
        end
        
        %low threshold
        cts=ctcounts(otherix);
        otherlow=lowthr(otherix);
        discard=isinf(otherlow)|isnan(otherlow)|otherlow==params.lowthr_clipval|cts'==0;
        cts(discard)=[];
        otherlow(discard)=[];
        if ~isempty(otherlow)
%             lowthr(i)=min(otherlow);  %most permissive
%             lowthr(i)=max(otherlow);  %most strict
            lowthr(i)=sum(cts/sum(cts).*otherlow(:)); %ct_count weigthed
        end
        
        %high threshold
        cts=ctcounts(otherix);
        otherhi=hithr(otherix);
        discard=isinf(otherhi)|isnan(otherhi)|otherhi==params.hithr_clipval|cts'==0;
        cts(discard)=[];
        otherhi(discard)=[];
        if ~isempty(otherhi)
%             hithr(i)=max(otherhi); %most permissive
%             hithr(i)=min(otherhi); %most strict
            hithr(i)=sum(cts/sum(cts).*otherhi(:)); %ct_count weigthed
        end
    end
    
    switch params.tailoption
        case 'low'
%             thisout=thisqc<=lowthr(i);
            hithr(i)=params.hithr_clipval; %remove from plots
        case 'high'
%             thisout=thisqc>=hithr(i);
            lowthr(i)=params.lowthr_clipval;
%         case 'both'
%             thisout=thisqc<=lowthr(i)|thisqc>=hithr(i);
    end
    thisout=thisqc<=lowthr(i)|thisqc>=hithr(i);
    
    outsub(i,thissub)=thisout;
    outliers(thissub)=outliers(thissub)|thisout;
end

result.params=params; 
result.vals=qcvals;
result.subix=tf;
result.outsub=outsub;
result.outliers=outliers;
result.lowthr=lowthr;
result.hithr=hithr;
result.mean=means;
result.median=medians;
result.std=stds;
result.mad0=mad0s;
result.mad1=mad1s;
result.geomean=geomeans;
result.geostd=geostds;