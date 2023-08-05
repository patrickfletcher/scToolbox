function [norm, hvg, mnn] = scran_fastMNN(datafile, cellsub, splitby, genes, normpars, hvgpars, mnnpars, options)
arguments
    datafile
    cellsub
    splitby {string,char,cellstr}
    genes
    normpars.gene_subset=false
    normpars.do_multibatch=true
    normpars.do_pooledsizefactors=false
    normpars.min_mean=0.1 %key param
    hvgpars.n_features=0 %return all passing other thresholds
    hvgpars.do_poissonvar=false
    hvgpars.do_densityweights=false %overrides fitTrendVar default of TRUE. important when HVGs are also high-abundance
    hvgpars.var_thr=0.0 %key param
    hvgpars.fdr_thr=1  % set to 1 to omit FDR threshold. overly conservative.
    hvgpars.min_mean_hvg=0.1 %key param
    mnnpars.k=20
    mnnpars.prop_k=[]
    mnnpars.d=50
    mnnpars.ndist=3
    mnnpars.merge_order=[]
    options.tmp_path="D:/tmp/tmp_scran_fastMNN/"
    options.tmp_fileroot="tmp"
%     options.resultfile="D:/tmp/tmp_scran_fastMNN_results.mat"
    options.verbose=false
end
% TODO - complete api:
% - control size factors
% - control HVG selection
%TODO: save the SCE object

% do_poissonvar = as.logical(hvgpars$do.poissonvar)
% do_topn = as.logical(hvgpars$do.topn)
% n_hvg = as.numeric(hvgpars$n.features)
% var_thr = as.numeric(hvgpars$var.thr)
% fdr_thr = as.numeric(hvgpars$fdr.thr)

if options.verbose
    disp("Running " + mfilename + "...")
end

if ~exist(options.tmp_path,"dir")
    mkdir(options.tmp_path)
end

datafile=cellstr(datafile);
cellsubsetfile='D:/tmp/tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)

cellsubsetfile=cellstr(cellsubsetfile);
splitby=cellstr(splitby);

parfile=fullfile(options.tmp_path, options.tmp_fileroot+"_pars.mat");
save(parfile,"datafile","cellsubsetfile","splitby","normpars","hvgpars","mnnpars")

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


norm=normpars;
sfsfile=fullfile(options.tmp_path, options.tmp_fileroot+"_sfs.txt");
norm.sizefactors=readmatrix(sfsfile);

hvg=hvgpars;
hvgfile=fullfile(options.tmp_path, options.tmp_fileroot+"_hvgs.txt");
hvg.id=string(readcell(hvgfile));
hvg.ix= getGeneIndices(hvg.id,genes.id);
hvg.name = genes.name(hvg.ix);

mnn=mnnpars;
mnnfile=fullfile(options.tmp_path, options.tmp_fileroot+"_mnn.csv");
rotfile=fullfile(options.tmp_path, options.tmp_fileroot+"_rot.csv");
infofile=fullfile(options.tmp_path, options.tmp_fileroot+"_minfo.csv");
mnn.coords=readmatrix(mnnfile);
mnn.coeff=readmatrix(rotfile);
mnn.merge_info = readtable(infofile);

clustfile=fullfile(options.tmp_path, options.tmp_fileroot+"_qclust.txt");
if normpars.do_pooledsizefactors
    norm.clust=string(readcell(clustfile));
end
toc