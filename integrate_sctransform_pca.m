function [res, result_tab] = integrate_sctransform_pca(datafile, cellsub, splitby, genesub, options)
arguments
    datafile
    cellsub
    splitby
    genesub=1
    options.n_pcs = 50
    options.do_compute = true;
    options.result_file = 'tmp_sct_pca.csv'
    options.verbose = false
    options.saveSO = 0
    options.saveSOtag = ""
end


%cellsub must be a table with two columns: id, keep, <splitby grouping variable>
cellsubfile='tmp_sct_cellsub.csv';
writetable(cellsub,cellsubfile)

genesubfile='tmp_sct_genesub.csv';
if ~isnumeric(genesub)
    writetable(genesub,genesubfile)
else
    genesubfile=genesub;
end

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\integrate_sctransform_pca.R";
command = char(strjoin([scriptfile, datafile, cellsubfile, options.n_pcs, ...
    options.result_file, splitby, options.saveSO, options.saveSOtag, genesubfile]," "));

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

result_tab=readtable(options.result_file);

res.npc=options.n_pcs;
res.method="SCT";
% res.coords=PCs;
res.coords=result_tab{:,2:end};

[path,name,ext]=fileparts(options.result_file);
res.loadings=readtable(fullfile(path,"feature_loadings_"+name+ext));
res.stdev=readtable(fullfile(path,"stdev_"+name+ext));