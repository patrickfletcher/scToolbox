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
    
    lowthr(i)=max(lowthr(i),params.lowthr_clipval); %clamp thresholds
    hithr(i)=min(hithr(i),params.hithr_clipval); 
    
    
    if ismember(thisctname,params.manual_low.names)
        lowthr(i)=params.manual_low.vals(params.manual_low.names==thisctname);
    end
    if ismember(thisctname,params.manual_high.names)
        hithr(i)=params.manual_high.vals(params.manual_high.names==thisctname);
    end
    
    %special value for "Unc": lowthr=min(others), hithr=max(others)
    % Assumes Unc is last. 
    % alt versions? 
    % - cell count weigthed thr
    % - min/max of kept values observed in any other type
    if ctnames(i)=="Unc"
        otherlow=lowthr(1:i-1);
        cts=ctcounts(1:i-1);
        cts(isnan(otherlow))=[];
        cts(otherlow==params.lowthr_clipval)=[];
%         otherlow(ctnames(1:i-1)=="T")=[];
        otherlow(isnan(otherlow))=[];
        otherlow(otherlow==params.lowthr_clipval)=[];
        if ~isempty(otherlow)
%             lowthr(i)=min(otherlow); 
            lowthr(i)=sum(cts/sum(cts).*otherlow(:));  
        end
        
        otherhi=hithr(1:i-1);
        cts=ctcounts(1:i-1);
        cts(isnan(hithr))=[];
        cts(otherhi==params.hithr_clipval)=[];
%         otherhi(ctnames(1:i-1)=="T")=[];
        otherhi(isnan(otherhi))=[];
        otherhi(otherhi==params.hithr_clipval)=[];
        if ~isempty(otherhi)
%             hithr(i)=max(otherhi); 
            hithr(i)=sum(cts/sum(cts).*otherhi(:)); 
        end
    end
    
    switch params.tailoption
        case 'low'
            thisout=thisqc<=lowthr(i);
            hithr(i)=nan; %remove from plots
        case 'high'
            thisout=thisqc>=hithr(i);
            lowthr(i)=nan;
        case 'both'
            thisout=thisqc<=lowthr(i)|thisqc>=hithr(i);
    end
    
    outsub(i,thissub)=thisout;
    outliers(thissub)=outliers(thissub)|thisout;
end

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