function T=compile_summaries(filepath, outputfile)

listing=dir(filepath);
listingnames={listing.name};
iscsv=contains(listingnames,'.csv');
csvfiles=listingnames(iscsv);

opts = detectImportOptions(fullfile(filepath,csvfiles{1}));
opts = setvartype(opts,'string');

T=table();
names=[];
for i=1:length(csvfiles)
    csvfile=csvfiles{i};
    [~,name]=fileparts(csvfile);
    names=[names;string(name)];
    T=[T;readtable(fullfile(filepath,csvfile),opts)];
end
N=table(names,'VariableNames',{'name'});
T=[N,T];

if exist('outputfile','var')
    writetable(T,fullfile(filepath,outputfile));
end