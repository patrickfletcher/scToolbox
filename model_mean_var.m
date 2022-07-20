function result =  model_mean_var(expr, options)
arguments
    expr
    options.min_mean=0.1
    options.gene_subset=[] %only consider this set of genes overall
    options.fit_subset=[] %use only these to fit the technical component
    options.block=[]
    options.weight_blocks=false
end
% modeled from SCRAN.
% - fit a trend curve (LOWESS) to mean-variance of log-expression data. Assuming
% most (enough?) genes are not variable, this fit should reflect the
% uninteresting (technical) part of the mean-variance.
% - subtract the trend: residuals are "biological component"
% - need method (test?) for: is this gene more variable than other genes of the same abundance?
% -- scran assumes null gene residuals are normally distributed about 0
% with spread as measured during fitting...

% they first fit: y = ax/(x^n + b)
% then do something with lowess..
%trend fitting: use "smoothdata" - lowess, loess, rlowess, rloess