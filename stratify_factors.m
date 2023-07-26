function [result, name_combos]=stratify_factors(factor1,factor2, opts)
arguments
    factor1
    factor2
    opts.separator="_"
    opts.min_count=1
    opts.removecats = true
end
%combine two categorical arrays into one stratified array

if ~iscategorical(factor1)
    factor1=categorical(factor1);
end
if ~iscategorical(factor2)
    factor2=categorical(factor2);
end

separator=string(opts.separator);

%remove cats first?
if opts.removecats
    factor1=removecats(factor1);
    factor2=removecats(factor2);
end

factor1Names=categories(factor1);
factor2Names=categories(factor2);
result=factor1;
result(:)=missing(); result=removecats(result);
for i=1:length(factor1Names)
    thisf1=factor1==factor1Names{i};
    for j=1:length(factor2Names)
        thisf2=factor2==factor2Names{j};
        result(thisf1&thisf2)=factor1Names{i}+separator+factor2Names{j};
    end
end

cats = categories(result);
C=countcats(result);
for i=1:length(cats)
    if C(i)<opts.min_count
        result(result==cats{i})=missing();
    end
end
if opts.removecats
    result = removecats(result);
end