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
outsub=false(size(tf));
outliers=false(size(qcvals));
for i=1:size(tf,1)
    thissub=tf(i,:);
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
            hithr(i)=prctile(thisqc,params.hiptile);
            
    end
    
    lowthr(i)=max(lowthr(i),params.lowthr_clipval); %clamp thresholds
    hithr(i)=min(hithr(i),params.hithr_clipval); 
    
    %special value for "Unc": lowthr=min(others), hithr=max(others)
    % Assumes Unc is last. 
    % alt versions? 
    % - cell count weigthed thr
    % - min/max of kept values observed in any other type
    if ctnames(i)=="Unc"
        otherlow=lowthr(1:i-1);
%         otherlow(ctnames(1:i-1)=="T")=[];
        otherlow(isnan(otherlow))=[];
        otherlow(otherlow==params.lowthr_clipval)=[];
        if ~isempty(otherlow)
            lowthr(i)=min(otherlow); 
        end
        otherhi=hithr(1:i-1);
%         otherhi(ctnames(1:i-1)=="T")=[];
        otherhi(isnan(otherhi))=[];
        otherhi(otherhi==params.hithr_clipval)=[];
        if ~isempty(otherhi)
            hithr(i)=max(otherhi); 
        end
    end
    
    for j=1:length(params.manual_low.names)
        lowthr(ctnames==params.manual_low.names(j))=params.manual_low.vals(j);
    end
    for j=1:length(params.manual_high.names)
        hithr(ctnames==params.manual_high.names(j))=params.manual_high.vals(j);
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