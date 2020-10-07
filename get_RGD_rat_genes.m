%get the latest RGD Rat Gene annotations and store as a table

rgdSaveFile='D:\scRNAseq\GENES_RAT.txt';
rgdSavePath=fileparts(rgdSaveFile);

% readRGDfromWeb
ftpobj = ftp('ftp.rgd.mcw.edu');
cd(ftpobj,'pub/data_release');
mget(ftpobj,'GENES_RAT.txt',rgdSavePath);

opts=detectImportOptions(rgdSaveFile,'TextType','string');
chrvar=opts.VariableNames(contains(opts.VariableNames,'CHROMOSOME'));
opts=setvartype(opts,chrvar,'string');
idvar=opts.VariableNames(contains(opts.VariableNames,'_ID'));
opts=setvartype(opts,idvar,'string');
RGD=readtable(rgdSaveFile,opts);

%subset of columns of interest
rgd=table();
rgd.NCBI_ID=RGD.NCBI_GENE_ID;
rgd.RGD_ID=RGD.GENE_RGD_ID;
rgd.ENS_ID=RGD.ENSEMBL_ID;
rgd.RGD_Symbol=RGD.SYMBOL;
rgd.name=RGD.NAME;
rgd.description=RGD.GENE_DESC;
rgd.gene_type=RGD.GENE_TYPE;
rgd.RGD_old_symbol=RGD.OLD_SYMBOL;
rgd.old_symbol=RGD.OLD_SYMBOL;
rgd.old_name=RGD.OLD_NAME;

save RGDtable.mat rgd