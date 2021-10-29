function [res, result_tab] = integrate_sctransform_pca(datafile, cellsub, splitby, options)
arguments
    datafile
    cellsub
    splitby
    options.n_pcs = 50
    options.do_compute = true;
    options.result_file = 'tmp_sct_pca.csv'
    options.verbose = false
    options.saveSO = 0
    options.saveSOtag = ""
end

cellsubfile='tmp_sct_cellsub.csv';

%cellsub must be a table with two columns: id, keep, <splitby grouping variable>
tic
writetable(cellsub,cellsubfile)
toc

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\integrate_sctransform_pca.R";
command = char(strjoin([scriptfile, datafile, cellsubfile, options.n_pcs, ...
    options.result_file, splitby, options.saveSO, options.saveSOtag]," "));

commandline=['"', Rpath, command '"'];

if options.do_compute %else, uses result_file
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


tic
result_tab=readtable(options.result_file);
toc
result_tab=renamevars(result_tab,"Var1","barcode");

res.npc=options.n_pcs;
res.method="SCT";
res.coords=result_tab{:,2:end};