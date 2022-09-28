function result = scran_fastMNN(datafile, cellsub, splitby, normpars, hvgpars, mnnpars, options)
arguments
    datafile
    cellsub
    splitby {string,char,cellstr}
    normpars.gene_subset=false
    normpars.do_multibatch=false
    normpars.do_pooledsizefactors=false
    normpars.min_mean=0.1
    hvgpars.n_features=3000
    hvgpars.do_poissonvar=true
    hvgpars.do_topn=false
    hvgpars.var_thr=0
    hvgpars.fdr_thr=1
    mnnpars.k=20
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

result.method="scran_fastMNN";
result.normpars=normpars;
result.hvgpars=hvgpars;
result.mnnpars=mnnpars;

datafile=cellstr(datafile);

cellsubsetfile='D:/tmp/tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

splitby=cellstr(splitby);
% mnnpars.merge_order
if isempty(mnnpars.merge_order)
%     mnnpars.merge_order=1:length(unique(cellsub.(string(splitby))));
%     mnnpars.merge_order=cellstr(unique(cellsub.(string(splitby))));
end

if ~exist(options.tmp_path,"dir")
    mkdir(options.tmp_path)
end

parfile=fullfile(options.tmp_path, options.tmp_fileroot+"_pars.mat");

save(parfile,"datafile","cellsubsetfile","splitby","normpars","hvgpars","mnnpars")

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
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

mnnfile=fullfile(options.tmp_path, options.tmp_fileroot+"_mnn.csv");
rotfile=fullfile(options.tmp_path, options.tmp_fileroot+"_rot.csv");
hvgfile=fullfile(options.tmp_path, options.tmp_fileroot+"_hvgs.txt");
result.coords=readmatrix(mnnfile);
result.coeff=readmatrix(rotfile);
result.hvgs=string(readcell(hvgfile));

% R_result=load(options.resultfile);
% result.coords=R_result.mnn;
% result.coeff=R_result.rot;
% result.hvgs=R_result.hvgs;
% result.sizefactors=R_result.sizefactors;