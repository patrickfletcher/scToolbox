function [res, result_tab] = sctransform_pca(datafile, cellsub, genesub, options)
arguments
    datafile
    cellsub
    genesub=1
    options.n_pcs = 50
    options.do_compute = true;
    options.use_so = false;
    options.result_file = 'tmp_sct_pca.csv'
    options.verbose = false
end
% if use_so==true, datafile is interpreted as the RDS file to load

%cellsub must be a table with two columns: id, keep
cellsubfile='tmp_sct_cellsub.csv';
writetable(cellsub,cellsubfile)

genesubfile='tmp_sct_genesub.csv';
if ~isnumeric(genesub)
    writetable(genesub,genesubfile)
else
    genesubfile=genesub;
end

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\sctransform_pca.R";
command = char(strjoin([scriptfile,datafile,cellsubfile,options.n_pcs,...
    options.result_file,options.use_so, genesubfile]," "));

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

result_tab=readtable(options.result_file);

res.npc=options.n_pcs;
res.method="SCT";
% res.coords=PCs;
res.coords=result_tab{:,2:end};

[path,name,ext]=fileparts(options.result_file);
res.loadings=readtable(fullfile(path,"feature_loadings_"+name+ext));
res.stdev=readtable(fullfile(path,"stdev_"+name+ext));