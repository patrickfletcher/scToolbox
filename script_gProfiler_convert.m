%script to convert gene sets from human to rat with gProfiler

%housekeeping genes (Tirosh16)
% HK = importdata('housekeepers_human.txt');
% save_file = 'housekeepers_human_rat.csv';
% genes = HK;

%cell cycle score
CC = importdata('regev_lab_cellcycle_human.txt');
save_file = 'regev_lab_cellcycle_human_rat.csv';
genes = CC;

gp=py.gprofiler.GProfiler();
ids=py.list(cellstr(genes')); %needs a row-vector cell array of chars
tic
gConvArgs=pyargs('organism','rnorvegicus','query',ids,'target_namespace','ENSG');
res=gp.convert(gConvArgs); %returns list of dicts
res=Core_py2matlab(res); 
res=[res{:}];
gEns=struct2table(res); %all IDs map, some gene symbols change.
toc

gTab = table();
gTab.human = gEns.incoming;
gTab.rat = gEns.name;
gTab.ens_id = gEns.converted;
gTab.description = gEns.description;
gTab.n_incoming = gEns.n_incoming;
gTab.n_converted = gEns.n_converted;
gTab.namespaces = gEns.namespaces;

writetable(gTab,save_file)
