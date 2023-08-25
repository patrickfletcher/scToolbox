function [norm, hvg] = scran_multibatch_hvg(datafile, cellsub, batchvar, genes, normpars, hvgpars, options)
arguments
    datafile
    cellsub
    batchvar {string,char,cellstr}
    genes
    normpars.gene_subset=false
    normpars.do_multibatch=true
    normpars.do_pooledsizefactors=false
    normpars.min_mean=0.1 %key param
    hvgpars.n_features=0
    hvgpars.do_poissonvar=false
    hvgpars.do_densityweights=false %overrides fitTrendVar default of TRUE. important when HVGs are also high-abundance
    hvgpars.var_thr=0.0 %key param
    hvgpars.fdr_thr=1 % set to 1 to omit FDR threshold. overly conservative.
    hvgpars.min_mean_hvg=0.1 %key param
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

if ~exist(options.tmp_path,"dir")
    mkdir(options.tmp_path)
end

datafile=cellstr(datafile);

cellsubsetfile='D:/tmp/tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

batchvar=cellstr(batchvar);

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

meanfile=fullfile(options.tmp_path, options.tmp_fileroot+"_mean.txt");
varfile=fullfile(options.tmp_path, options.tmp_fileroot+"_var.txt");
biofile=fullfile(options.tmp_path, options.tmp_fileroot+"_bio.txt");
techfile=fullfile(options.tmp_path, options.tmp_fileroot+"_tech.txt");
hvg.mean=readmatrix(meanfile);
hvg.var=readmatrix(varfile);
hvg.bio=readmatrix(biofile);
hvg.tech=readmatrix(techfile);

if normpars.do_pooledsizefactors
    clustfile=fullfile(options.tmp_path, options.tmp_fileroot+"_qclust.txt");
    norm.clust=string(readcell(clustfile));
end
toc