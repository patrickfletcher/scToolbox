library(R.matlab)

args <- commandArgs(TRUE)

parsmatfile <- args[1]
resultsmatfile <- args[2]

f <- readMat(parsmatfile)
#strings/chars should be cellstr

pars=f$pars[,,1]

print(pars)

res1=unlist(pars$par1)
res2=pars$par2

writeMat(resultsmatfile, result1=res1, result2=res2)