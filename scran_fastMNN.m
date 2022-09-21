function result = scran_fastMNN(datafile, cellsub, splitby, normpars, mnnpars, options)
arguments
    datafile
    cellsub
    splitby {string,char,cellstr}
    normpars.gene_subset=false
    normpars.n_features=3000
    normpars.do_multibatch=false
    normpars.do_pooledsizefactors=false
    normpars.min_mean=0.1
    mnnpars.k=20
    mnnpars.d=50
    mnnpars.ndist=3
    mnnpars.merge_order=[]
    options.parfile="D:/tmp/tmp_scran_fastMNN_pars.mat"
    options.resultfile="D:/tmp/tmp_scran_fastMNN_results.mat"
    options.verbose=false
end
% TODO - complete api:
% - control size factors
% - control HVG selection
%TODO: save the SCE object

if options.verbose
    disp("Running " + mfilename + "...")
end

result.method="scran_fastMNN";
result.normpars=normpars;
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

save(options.parfile,"datafile","cellsubsetfile","splitby","normpars","mnnpars")

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
result.coords=R_result.mnn;
result.coeff=R_result.rot;
result.hvgs=R_result.hvgs;
result.sizefactors=R_result.sizefactors;