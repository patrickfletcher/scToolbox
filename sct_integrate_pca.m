function res = sct_integrate_pca(datafile, cellinfo, splitby, sctpars, intpars, outopts, options)
arguments
    datafile
    cellinfo
    splitby
    sctpars.vars_to_regress='NULL'
    sctpars.variable_features_n=3000
    sctpars.batch_var='NULL'
    sctpars.n_genes=2000
    sctpars.n_cells=5000
    sctpars.bin_size=500
    sctpars.min_cells=5
    sctpars.vst_method='glmpoisson'
    sctpars.vst_flavor='NULL'
    intpars.n_anchor_features=2000
    intpars.n_pcs = 50
    intpars.reduction='rpca'
    intpars.k_anchor=5
    intpars.k_filter=200
    intpars.k_score=30
    intpars.max_features=200
    intpars.k_weight=100
    intpars.refs='NULL'
    intpars.sample_tree='NULL'
    outopts.out_pcs = 50
    outopts.save_sctlist = 0
    outopts.save_so = 0
    outopts.so_tag = ""
    options.parfile="tmp_sct_integrate_pca_pars.mat"
    options.resultfile="tmp_sct_integrate_pca_results.mat"
    options.verbose=true
end

% vars.to.regress = unlist(sctpars$vars.to.regress)
% variable.features.n = sctpars$variable.features.n
% batch.var = unlist(sctpars$batch.var)
% n.genes = sctpars$n.genes  # Number of genes to use when estimating parameters
% n.cells = sctpars$n.cells  # Number of cells to use when estimating parameters
% bin.size = sctpars$bin.size
% min.cells = sctpars$min.cells
% vst.method = unlist(sctpars$vst.method)
% vst.flavor = unlist(sctpars$vst.flavor) 

% nfeatures = intpars$n.anchor.features
% n.pcs = intpars$n.pcs
% reduction = unlist(intpars$reduction)
% k.anchor = intpars$k.anchor
% k.filter = intpars$k.filter
% k.score = intpars$k.score
% max.features = intpars$max.features
% k.weight = intpars$k.weight
% refs = intpars$refs
% sample.tree = intpars$sample.tree

% out.pcs = outopts$out.pcs
% save.so = outopts$save.so
% so.tag = outopts$so.tag

%make sure all char/string options are cellstr
datafile=cellstr(datafile);
splitby=cellstr(splitby);
sctpars.vars_to_regress=cellstr(sctpars.vars_to_regress);
sctpars.batch_var=cellstr(sctpars.batch_var);
sctpars.vst_method=cellstr(sctpars.vst_method);
sctpars.vst_flavor=cellstr(sctpars.vst_flavor);
intpars.reduction=cellstr(intpars.reduction);
if ~isnumeric(intpars.refs)
    intpars.refs=cellstr(intpars.refs);
end
if ~isnumeric(intpars.sample_tree)
    intpars.sample_tree=cellstr(intpars.sample_tree);
end
outopts.so_tag=cellstr(outopts.so_tag);

%csv file to handle the cellinfo table...
cellinfofile='tmp_sct_cellsub.csv';
writetable(cellinfo, cellinfofile)
cellinfofile=cellstr(cellinfofile);

save(options.parfile,"datafile","cellinfofile","splitby","sctpars","intpars","outopts")

tmpresultfile=options.resultfile+".tmp";
if exist(tmpresultfile,"file")
    delete(tmpresultfile)
end

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = mfilename('fullpath')+".R";
command = char(strjoin([scriptfile,options.parfile,options.resultfile]," "));
commandline=['"', Rpath, command '"'];

tic
[status, cmdout]=system(commandline);

if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end
disp(mfilename+": ")
toc

res.method="sct_integrate_pca";
res.sctpars=sctpars;
res.intpars=intpars;
res.npc=outopts.out_pcs;

R_result=load(options.resultfile);

res.coords=R_result.pca;
res.features=R_result.features;
res.coeff=R_result.coeff;
res.pc_std=R_result.pc_std;
