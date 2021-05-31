function genetab = get_biomart_genes(go_list, go_file, gene_file)
arguments
    go_list
    go_file = 'tmp_go_list.csv'
    gene_file = 'tmp_gene_list.csv'
end

writecell(go_list,go_file)

Rpath = "C:\Users\fletcherpa\Documents\R\R-4.0.5\bin\Rscript --vanilla ";
[status, cmdout]=system(Rpath + strjoin(["get_go_genes.R",go_file,gene_file]," "));
if status~=0
    error(cmdout)
end
genetab=readtable(gene_file);
genetab=sortrows(genetab);