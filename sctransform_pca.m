function [res, result_tab] = sctransform_pca(datafile, cellsub, options)
arguments
    datafile
    cellsub
    options.n_pcs = 50
    options.do_compute = true;
    options.use_so = false;
    options.result_file = 'tmp_sct_pca.csv'
    options.verbose = false
end
% if use_so==true, datafile is interpreted as the RDS file to load

cellsubfile='tmp_sct_cellsub.csv';

%cellsub must be a table with two columns: id, keep
tic
writetable(cellsub,cellsubfile)
toc

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\sctransform_pca.R";
command = char(strjoin([scriptfile,datafile,cellsubfile,options.n_pcs,...
    options.result_file,options.use_so]," "));

commandline=['"', Rpath, command '"'];

if options.do_compute
    tic
    [status, cmdout]=system(commandline);
    toc
    if status~=0
        error(cmdout)
    end
    if options.verbose
        disp(cmdout)
    end
end
%else, uses result_file

tic
result_tab=readtable(options.result_file);
toc
result_tab=renamevars(result_tab,"Var1","barcode");


res.npc=options.n_pcs;
res.method="SCT";
res.coords=result_tab{:,2:end};