function S = GOSemSim(term_ids, ontology, options)
arguments
    term_ids
    ontology %BP, CC, MF
    options.simmethod="Wang"
end

tmpfileroot = strjoin(["tmp_GOSemSim",ontology,options.simmethod],"_");
subfile = tmpfileroot+"_GOids.csv";
writecell(term_ids,subfile)

outfile=tmpfileroot+"_matrix.csv";

Sids={};
if exist(outfile,'file')
    Stab=readtable(outfile,'ReadVariableNames',false);
    Sids=Stab.Var1;
end
if ~isequal(Sids,term_ids(:))
    Rpath = "C:\Users\fletcherpa\Documents\R\R-4.1.0\bin\Rscript --vanilla ";
    script = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\GOSemSim.R";
    [status, cmdout]=system(Rpath + strjoin([script,subfile,ontology,options.simmethod,outfile]," "));
    if status~=0
        error(cmdout)
    end
    Stab=readtable(outfile,'ReadVariableNames',false);
end
S=Stab{:,2:end};