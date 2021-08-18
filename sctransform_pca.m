function [res, result_tab] = sctransform_pca(datafile, cellsub, options)
arguments
    datafile
    cellsub
    options.n_pcs = 50
    options.result_file = 'tmp_sct_pca.csv'
    options.verbose = false
end

cellsubfile='tmp_sct_cellsub.csv';

%cellsub must be a table with two columns: id, keep
writetable(cellsub,cellsubfile)

Rpath = "C:\Users\fletcherpa\Documents\R\R-4.1.0\bin\Rscript --vanilla ";
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\run_SCTransform_PCA.R";

[status, cmdout]=system(Rpath + strjoin([scriptfile,datafile,cellsubfile,options.n_pcs,options.result_file]," "));
if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end

result_tab=readtable(options.result_file);
result_tab=renamevars(result_tab,"Var1","barcode");


res.npc=options.n_pcs;
res.method="SCT";
res.coords=result_tab{:,2:end};