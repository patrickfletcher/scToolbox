function tfidf = quickMarkers(countfile, clusters, options)
arguments
    countfile
    clusters
    options.N=10
    options.FDR=0.05
    options.expressCut=0.9
    options.result_file='tmp_quickMarkers_result.csv'
    options.verbose=true
end

%cellsub must be a table with two columns: id, clusterID
clustfile='tmp_quickMarkers_clusters.csv';
writetable(clusters,clustfile)

Rpath = "C:\Users\fletcherpa\Documents\R\R-4.1.0\bin\Rscript --vanilla ";
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\quickMarkers.R";

[status, cmdout]=system(Rpath + strjoin([scriptfile,countfile,clustfile,...
    options.N,options.FDR,options.expressCut,options.result_file]," "));
if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end

tfidf=readtable(options.result_file);

