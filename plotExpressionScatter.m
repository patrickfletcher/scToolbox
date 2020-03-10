function [ax, hs, hc]=plotExpressionScatter(genelist,X,tcounts,genes,doRelative,figID,subplotdims,sp_params,docolorbar)
%wrapper for plotScatter that unpacks desired genes from tcounts/genes
% e.g.:
%  plotExpressionScatter("GH1", TSNE, tcounts, genes)

%defaults for this function
if ~exist('doRelative','var')||isempty(doRelative)
    doRelative=false;
end

%punt on defaults for plotScatter (let it handle those)
if ~exist('sp_params','var')||isempty(sp_params)
    sp_params=[];
end
if ~exist('subplotdims','var')||isempty(subplotdims)
    subplotdims=[];
end
if ~exist('figID','var')||isempty(figID)
    figID=[];
end
if ~exist('docolorbar','var')||isempty(docolorbar)
    docolorbar=[];
end

gix=getGeneIndices(genelist,genes.name);
C = tcounts(gix,:);
if doRelative
    C = C - genes.thr(gix);
end
[ax, hs, hc]=plotScatter(X,"value",genelist,C,figID,subplotdims,sp_params,docolorbar);