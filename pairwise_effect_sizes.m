function result = pairwise_effect_sizes(expr, groups1, groups2, groupopts, pars)
arguments
    expr
    groups1
    groups2
    groupopts
    pars
end
% compute effect sizes for pairwise combinations of groups1 and groups2 
%
%expects table of expression summarized by group: expr, prct
%alternate use case could be: pass in ncounts/tcounts, grouping