function result = sct_Harmony(datafile, cellsub, splitby, sctpars, harmonypars, options)
arguments
    datafile
    cellsub
    splitby {string,char,cellstr}
    sctpars.vars_to_regress{string,char,cellstr}="none"
    sctpars.n_features=3000
    sctpars.n_anchor_features=2000
    sctpars.n_pcs=50
    harmonypars.maxit=10
    harmonypars.saveit=false
    options.parfile="tmp_sct_Harmony_pars.mat"
    options.resultfile="tmp_sct_Harmony_results.mat"
    options.verbose=false
end
%would like to be able to just move forward one iteration at a time

result.method="sct_Harmony";
result.sctpars=sctpars;
result.harmonypars=harmonypars;

datafile=cellstr(datafile);

cellsubsetfile='tmp_sct_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

splitby=cellstr(splitby);

sctpars.vars_to_regress=cellstr(sctpars.vars_to_regress);

save(options.parfile,"datafile","cellsubsetfile","splitby","sctpars","harmonypars")

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
result.coords=R_result.harmony;