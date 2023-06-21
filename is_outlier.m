function [outliers, result] = is_outlier(data, type, group, options)
arguments
    data
    type = "both"
    group = [] %support batch-wise computations
    options.nmads=3
    options.do_log=false
    options.madconst=1.4826 %gives std scale
    options.min_diff=0
    options.lower_min=-inf;
    options.lower_max=inf;
    options.upper_min=-inf;
    options.upper_max=inf;
    options.clip_to_data=true
end
%compute median +/- k*mad threshold, identify outliers. 
% support vector of columns for data, with matching vectors of parameters

%simplify interface by not doing global clipping here. It's easy to do
%externally

%TODO: share missingness (for log?), share medians
%TODO: support providing subset of data to work on?

[nobs,nvars]=size(data);
if nobs<nvars
    error('data has fewer rows than columns: expects observations x feature in data')
end
if istable(data)
    varnames=data.Properties.VariableNames;
    data=data{:,:}; %no need for table anymore
else
    varnames="v"+string(1:nvars);
end

type=cellstr(type);
do_log=options.do_log;
nmads=options.nmads;
min_diff=options.min_diff;
lower_min=options.lower_min;
lower_max=options.lower_max;
upper_min=options.upper_min;
upper_max=options.upper_max;

if isempty(group)
    group=ones(nobs,1);
else
    if length(group)~=nobs
        error('group variable has different length than data observations')
    end
end
[g,gN]=grp2idx(group);

%expand parameters into vectors if needed (same thing for all vars)
if length(type)==1
    type=repmat(type,1,nvars);
end
if length(do_log)==1
    do_log=repmat(do_log,1,nvars);
end
if length(nmads)==1
    nmads=repmat(nmads,1,nvars);
end
if length(min_diff)==1
    min_diff=repmat(min_diff,1,nvars);
end

if length(lower_min)==1
    lower_min=repmat(lower_min,1,nvars);
end
if length(lower_max)==1
    lower_max=repmat(lower_max,1,nvars);
end
if length(upper_min)==1
    upper_min=repmat(upper_min,1,nvars);
end
if length(upper_max)==1
    upper_max=repmat(upper_max,1,nvars);
end

%adjust clip values in case of log
% lower_min(do_log&lower_min==-inf)=0;
% upper_max(do_log&upper_max==-inf)=0;

%compute the stats
result=struct('name',{},'type',{},'nmads',{},'do_log',{}, ...
    'lower_min',{},'lower_max',{},'upper_min',{},'upper_max',{},...
    'stats',{},'outlier',{});
outliers = false(size(group));
for i=1:nvars
    thisdata=data(:,i);
    if do_log(i)
        thisdata=log(thisdata);
    end
    meds = splitapply(@median,thisdata,g);
    mads = splitapply(@(x)mad(x,1),thisdata,g); %median abs dev
    mads = mads * options.madconst;
    
    %compute the thresholds
    diffval = max(nmads(i)*mads, min_diff(i));
    lower = meds - diffval;
    upper = meds + diffval;
    
    % remove unused constraints for upper/lower
    switch type{i}
        case "lower"
            upper = inf(size(upper));
        case "upper"
            lower = -inf(size(lower));
    end

    %clip to data range if desired
    for j=1:length(gN)
        if options.clip_to_data
            lower(j) = max(lower(j), min(thisdata(g==j)));
            upper(j) = min(upper(j), max(thisdata(g==j)));
        end
    end

    lower_noclip=lower;
    upper_noclip=upper;

    %apply global clip values
    if do_log(i) 
        lower = min(lower, log(lower_max(i)));
        lower = max(lower, log(lower_min(i)));
        upper = min(upper, log(upper_max(i)));
        upper = max(upper, log(upper_min(i)));
    else
        lower = min(lower, lower_max(i));
        lower = max(lower, lower_min(i));
        upper = min(upper, upper_max(i));
        upper = max(upper, upper_min(i));
    end

    this_outliers = false(size(group));
    for j=1:length(gN)
        this_outliers(g==j)=thisdata(g==j)<lower(j) | thisdata(g==j)>upper(j);
    end

    if do_log(i)
        meds=exp(meds);
%         mads=2.^mads;
%         diffval=2.^diffval;
        lower=exp(lower);
        upper=exp(upper);
        lower_noclip=exp(lower_noclip);
        upper_noclip=exp(upper_noclip);
    end

    %pack the results.
    res.name=string(varnames{i});
    res.type=string(type{i});
    res.nmads=nmads(i);
    res.do_log=do_log(i);
    res.lower_min=lower_min(i);
    res.lower_max=lower_max(i);
    res.upper_min=upper_min(i);
    res.upper_max=upper_max(i);
    stats=table();
    stats.groupnames=gN;
    stats.med=meds;
%     stats.mad=mads;
%     stats.diffval=diffval;
    stats.lower=lower;
    stats.upper=upper;
    stats.lower_noclip=lower_noclip;
    stats.upper_noclip=upper_noclip;
    res.stats=stats;
    res.outlier=this_outliers;
    result(end+1)=res;

    outliers=outliers | this_outliers;
end

