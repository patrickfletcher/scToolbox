function result = scran_hvgs(datafile, cellsub, genes, batch, normpars, hvgpars, options)
arguments
    datafile
    cellsub
    genes
    batch {string,char,cellstr}
    normpars.do_multibatch=true
    normpars.do_pooledsizefactors=true
    normpars.min_mean=0.1
    hvgpars.n_features=3000
    options.parfile="D:\tmp\tmp_scran_hvg_pars.mat"
    options.resultfile="D:\tmp\tmp_scran_hvg_results.mat"
    options.verbose=false
end
% TODO - complete api:
% - control size factors
% - control HVG selection
%TODO: save the SCE object

if options.verbose
    disp("Running " + mfilename + "...")
end

result.method="scran_hvgs";
result.normpars=normpars;
result.hvgpars=hvgpars;

datafile=cellstr(datafile);

cellsubsetfile='D:\tmp\tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

batch=cellstr(batch);

save(options.parfile,"datafile","cellsubsetfile","batch","normpars","hvgpars")

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

R_result=load(options.resultfile);
result.name=R_result.hvgs;
result.ix = getGeneIndices(R_result.hvgs, genes.name);
result.sizefactors=R_result.sizefactors;