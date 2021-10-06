function eulerr(sets, setnames, pdf_file, options)
arguments
    sets
    setnames
    pdf_file
    options.shape = "circle"
    options.cols = []
    options.alpha = 1
    options.fsz_counts = 8
    options.fsz_labels = 8
    options.width = 0
    options.height = 0
    options.rng_seed = 42
    options.verbose=false
end
%set fsz to 0 to not display the counts or labels

if ~iscell(sets) && isnumeric(sets)
    counts=sets(:);
    int_names=setnames(:);
else
    setnames=string(setnames);
    [excl_int,~,int_names]=get_exclusive_intersections(sets, setnames);
    int_names=string(int_names);
    int_names=strrep(int_names,'_','&');
    counts=cellfun(@length,excl_int);
end

T=table();
T.names=int_names(:);
T.counts=counts(:);

if ~isempty(options.cols)
    hexcols=rgb2hex(options.cols);
    hexcols=string(hexcols);
    T.cols=hexcols(:);
end

T(counts==0,:)=[];

setsfile="tmp_eulerr_data.csv";
writetable(T, setsfile);

Rpath = "C:\Users\fletcherpa\Documents\R\R-4.1.0\bin\Rscript --vanilla ";
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\eulerr.R";

pdf_file=fullfile(pwd,pdf_file);
[status, cmdout]=system(Rpath + strjoin([scriptfile,setsfile,options.shape,options.alpha,...
    options.width,options.height,options.fsz_counts,options.fsz_labels,options.rng_seed,pdf_file]," "));
if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end