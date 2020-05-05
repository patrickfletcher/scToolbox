function result=stratify_factors(factor1,factor2)
%combine two categorical arrays into one stratified array

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
        result(thisf1&thisf2)=factor1Names{i}+"_"+factor2Names{j};
    end
end