function result = test_R_matlab(pars, options)
arguments
    pars.par1 {string,char,cellstr}="test"
    pars.par2=[15,1.5]
    options.parfile="testpars.mat"
    options.resultfile="testresults.mat"
    options.verbose=true
end

pars.par1=cellstr(pars.par1);

save(options.parfile,"pars")

tmpresultfile=options.resultfile+".tmp";
if exist(tmpresultfile,"file")
    delete(tmpresultfile)
end

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = mfilename('fullpath')+".R";
command = char(strjoin([scriptfile,options.parfile,options.resultfile]," "));
commandline=['"', Rpath, command '"'];

[status, cmdout]=system(commandline);

if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end

result=load(options.resultfile);