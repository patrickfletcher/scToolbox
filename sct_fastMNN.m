function result = sct_fastMNN(datafile, cellsub, splitby, sctpars, mnnpars, options)
arguments
    datafile
    cellsub
    splitby {string,char,cellstr}
    sctpars.vars_to_regress{string,char,cellstr}="NULL"
    sctpars.n_features=3000
    sctpars.n_anchor_features=2000
    mnnpars.k=20
    mnnpars.d=50
    mnnpars.ndist=3
    mnnpars.merge_order=[]
    options.parfile="tmp_sct_fastMNN_pars.mat"
    options.resultfile="tmp_sct_fastMNN_results.mat"
    options.verbose=false
end

result.method="SCT_fastMNN";
result.sctpars=sctpars;
result.mnnpars=mnnpars;

datafile=cellstr(datafile);

cellsubsetfile='tmp_sct_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

splitby=cellstr(splitby);

sctpars.vars_to_regress=cellstr(sctpars.vars_to_regress);

% mnnpars.merge_order
if isempty(mnnpars.merge_order)
%     mnnpars.merge_order=1:length(unique(cellsub.(string(splitby))));
%     mnnpars.merge_order=cellstr(unique(cellsub.(string(splitby))));
end


save(options.parfile,"datafile","cellsubsetfile","splitby","sctpars","mnnpars")

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