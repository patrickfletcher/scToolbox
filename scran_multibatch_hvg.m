function [norm, hvg] = scran_multibatch_hvg(datafile, cellsub, batchvar, genes, normpars, hvgpars, options)
arguments
    datafile
    cellsub
    batchvar {string,char,cellstr}
    genes
    normpars.gene_subset=false
    normpars.do_multibatch=true
    normpars.do_pooledsizefactors=false
    normpars.min_mean=0.1
    hvgpars.n_features=[]
    hvgpars.do_poissonvar=false
    hvgpars.do_densityweights=false
    hvgpars.var_thr=0.001
    hvgpars.fdr_thr=1
    options.tmp_path=[]
    options.tmp_fileroot="tmp"
    options.verbose=false
end
%TODO: support user clusters for pooled sfs
% gene subset? aka row_subset in scran?

% do_poissonvar = as.logical(hvgpars$do.poissonvar)
% do_topn = as.logical(hvgpars$do.topn)
% n_hvg = as.numeric(hvgpars$n.features)
% var_thr = as.numeric(hvgpars$var.thr)
% fdr_thr = as.numeric(hvgpars$fdr.thr)

if options.verbose
    disp("Running " + mfilename + "...")
end

if isempty(options.tmp_path)
    options.tmp_path="D:/tmp/tmp_"+mfilename+"/";
end


hvgpars.do_topn=~isempty(hvgpars.n_features);

datafile=cellstr(datafile);

cellsubsetfile='D:/tmp/tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

batchvar=cellstr(batchvar);

if ~exist(options.tmp_path,"dir")
    mkdir(options.tmp_path)
end

parfile=fullfile(options.tmp_path, options.tmp_fileroot+"_pars.mat");

save(parfile,"datafile","cellsubsetfile","batchvar","normpars","hvgpars")

Rpath = [getExternalPath("R"), filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = mfilename('fullpath')+".R";
command = char(strjoin([scriptfile,options.tmp_path, options.tmp_fileroot]," "));
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

%TODO: link actual methods - pooledSFS vs libsize (pre multibatch)
normpars.method="scran_multibatchnorm";
hvgpars.method="scran_hvgs";

norm=normpars;
hvg=hvgpars;

%TODO: return regular libsizefactors
sfsfile=fullfile(options.tmp_path, options.tmp_fileroot+"_sfs.txt");
norm.sizefactors=readmatrix(sfsfile);
% norm.libsizefactors

%TODO: pass gene table in to get name/ix
hvgfile=fullfile(options.tmp_path, options.tmp_fileroot+"_hvgs.txt");
hvg.id=string(readcell(hvgfile));
hvg.ix= getGeneIndices(hvg.id,genes.id);
hvg.name = genes.name(hvg.ix);

if normpars.do_pooledsizefactors
    clustfile=fullfile(options.tmp_path, options.tmp_fileroot+"_qclust.txt");
    norm.clust=string(readcell(clustfile));
end