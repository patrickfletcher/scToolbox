function result = scran_singleR(datafile, refinfo, cellinfo, normpars, options)
arguments
    datafile
    refinfo
    cellinfo
    normpars.n_features=3000
    normpars.do_multibatch=true
    normpars.do_pooledsizefactors=true
    normpars.min_mean=0.1
    options.parfile="D:\tmp\tmp_scran_fastMNN_pars.mat"
    options.resultfile="D:\tmp\tmp_scran_fastMNN_results.mat"
    options.verbose=false
end