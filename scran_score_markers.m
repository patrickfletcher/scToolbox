function result = scran_score_markers(datafile, cellsub, groupvar, batchvar, normpars, scorepars, options)
arguments
    datafile
    cellsub
    groupvar {string,char,cellstr}
    batchvar {string,char,cellstr} = "NULL"
    normpars.do_multibatch_norm=false
    normpars.do_pooledsizefactors=false
    normpars.min_mean=0.1
    scorepars.pairings="NULL"  %self
    scorepars.pairings2="NULL"  %other
    scorepars.lfc=0
    scorepars.full_stats=false
    scorepars.min_cells=1
    options.parfile="D:\tmp\tmp_scran_score_markers.mat"
    options.resultfile="D:\tmp\tmp_scran_score_markers.mat"
    options.verbose=false
end

if options.verbose
    disp("Running " + mfilename + "...")
end

datafile=cellstr(datafile);

cellsubsetfile='D:\tmp\tmp_cellsub.csv';
writetable(cellsub,cellsubsetfile)
cellsubsetfile=cellstr(cellsubsetfile);

groupvar=cellstr(groupvar);
batchvar=cellstr(batchvar);

scorepars.pairings=cellstr(scorepars.pairings);
scorepars.pairings2=cellstr(scorepars.pairings2);

save(options.parfile,"datafile","cellsubsetfile","groupvar","batchvar","normpars","scorepars")

tmpresultfile=options.resultfile+".tmp";
if exist(tmpresultfile,"file")
    delete(tmpresultfile)
end

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = mfilename('fullpath')+".R";
command = char(strjoin([scriptfile,options.parfile,options.resultfile]," "));
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
toc

result.method="score_markers";
normfields=fieldnames(normpars);
for i=1:length(normfields)
    result.(normfields{i})=normpars.(normfields{i});
end
scorefields=fieldnames(scorepars);
for i=1:length(scorefields)
    result.(scorefields{i})=scorepars.(scorefields{i});
end

R_result=load(options.resultfile);
gnames=scorepars.pairings;
% gnames=cellstr(R_result.groupnames);
genenames=R_result.genenames;
genetab=table; genetab.gene=genenames(:);

result.mean=genetab; 
result.mean.Properties.RowNames=genenames(:);
means=structfun(@(x)double(x(:)),R_result.averages,'UniformOutput',false);
means=struct2table(means);

result.prop=genetab; 
result.prop.Properties.RowNames=genenames(:);
props=structfun(@(x)double(x(:)),R_result.props,'UniformOutput',false);
props=struct2table(props);

allnames=fieldnames(R_result.averages);
for i=1:length(allnames)
    result.mean.(allnames{i})=means.(allnames{i});
    result.prop.(allnames{i})=props.(allnames{i});
end

for i=1:length(gnames)
    thisG=R_result.(gnames{i});
    thisG=structfun(@(x)double(x(:)),thisG,'UniformOutput',false);
    thisG=struct2table(thisG);
    thisG.Properties.RowNames=genenames(:);
    result.(gnames{i})=[genetab,thisG];
end