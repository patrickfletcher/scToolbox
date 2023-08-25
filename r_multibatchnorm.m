function norm = r_multibatchnorm(datafile, cellsub, batchvar, normpars, options)
arguments
    datafile
    cellsub
    batchvar {string,char,cellstr}
    normpars.gene_subset=false
    normpars.min_mean=0.1 %key param - can be a list
    options.tmp_path="D:/tmp/tmp_"+mfilename+"/"
    options.tmp_fileroot="tmp"
    options.verbose=false
end
%TODO: support user clusters for pooled sfs
% gene subset? aka row_subset in scran?

if options.verbose
    disp("Running " + mfilename + "...")
end

if ~exist(options.tmp_path,"dir")
    mkdir(options.tmp_path)
end

datafile=cellstr(datafile);
cellsubsetfile=options.tmp_path + options.tmp_fileroot + "cellsub.csv";
writetable(cellsub,cellsubsetfile)

cellsubsetfile=cellstr(cellsubsetfile);
batchvar=cellstr(batchvar);

parfile=fullfile(options.tmp_path, options.tmp_fileroot+"_pars.mat");
save(parfile,"datafile","cellsubsetfile","batchvar","normpars")

Rpath = [getExternalPath("R"), filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = mfilename('fullpath')+".R";
command = char(strjoin([scriptfile,options.tmp_path, options.tmp_fileroot]," "));
commandline=['"', Rpath, command '"'];

tic
[status, cmdout]=system(commandline);

if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end
disp(mfilename+": ")

%TODO: link actual methods - pooledSFS vs libsize (pre multibatch)
normpars.method="multibatchnorm";

norm=normpars;

%TODO: return regular libsizefactors
sfsfile=fullfile(options.tmp_path, options.tmp_fileroot+"_sfs.txt");
norm.sizefactors=readmatrix(sfsfile);
toc