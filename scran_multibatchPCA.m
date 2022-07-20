function result = scran_multibatchPCA(datafile, cellsub, batchvar, normpars, pcapars, options)
arguments
    datafile
    cellsub
    batchvar {string,char,cellstr}
    normpars.n_features=3000
    normpars.do_multibatch_norm=true
    normpars.do_pooledsizefactors=true
    normpars.min_mean=0.1
    pcapars.d=50
    pcapars.w='NULL'
    options.parfile="D:\tmp\tmp_scran_multibatchPCA_pars.mat"
    options.resultfile="D:\tmp\tmp_scran_multibatchPCA_results.mat"
    options.verbose=false
end

if options.verbose
    disp("Running " + mfilename + "...")
end

result.method=mfilename;
result.normpars=normpars;
result.pcapars=pcapars;

datafile=cellstr(datafile);

cellsubsetfile='D:\tmp\tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

batchvar=cellstr(batchvar);

save(options.parfile,"datafile","cellsubsetfile","batchvar","normpars","pcapars")

% Rscript dispatch command
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

R_result=load(options.resultfile);

result.coords=R_result.pca;
result.rotation=R_result.pca_extra.rotation;
result.centers=R_result.pca_extra.centers;
result.hvgs=R_result.hvgs;
result.sizefactors=R_result.sizefactors;