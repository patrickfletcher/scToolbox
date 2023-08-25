function [norm, hvg, mnn] = scran_fastMNN(datafile, csubtab, splitby, genes, normpars, hvgpars, mnnpars, options)
arguments
    datafile
    csubtab
    splitby {string,char,cellstr}
    genes
    normpars.gene_subset=false
    normpars.min_mean=0.1 
    normpars.do_multibatch=true
    normpars.do_pooledsizefactors=false
    hvgpars.min_mean_hvg=0.1
    hvgpars.n_features=0     
    hvgpars.do_poissonvar=false
    hvgpars.do_densityweights=false %overrides fitTrendVar default of TRUE. important when HVGs are also high-abundance
    hvgpars.var_thr=0.0 
    hvgpars.fdr_thr=1
    mnnpars.merge_order=[]
    mnnpars.d=50
    mnnpars.k=20
    mnnpars.prop_k=[]
    mnnpars.ndist=3
    options.tmp_path="D:/tmp/tmp_"+mfilename+"/"
    options.tmp_fileroot="tmp"
    options.verbose=false
end
%TODO: save the SCE object to RDS?

if options.verbose
    disp("Running " + mfilename + "...")
end
if ~exist(options.tmp_path,"dir")
    mkdir(options.tmp_path)
end

datafile=cellstr(datafile);
cellsubsetfile=options.tmp_path + options.tmp_fileroot + "cellsub.csv";
writetable(csubtab,cellsubsetfile)

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