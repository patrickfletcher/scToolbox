function result=stratify_factors(factor1,factor2, separator)
arguments
    factor1
    factor2
    separator="_"
end
%combine two categorical arrays into one stratified array

if ~iscategorical(factor1)
    factor1=categorical(factor1);
end
if ~iscategorical(factor2)
    factor2=categorical(factor2);
end

separator=string(separator);

%remove cats first?
factor1=removecats(factor1);
factor2=removecats(factor2);

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