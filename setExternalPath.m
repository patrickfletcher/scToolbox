function setExternalPath(programs, path)
arguments
    programs = ["python";"R"]
    path = []
end
% create/modify a .env file to store paths to python and R (others?)
% it will be a CSV, for easy use as a table

sctoolboxpath=fileparts(mfilename('fullpath'));
pathsFile = fullfile(sctoolboxpath,'extPaths.csv');

%make the file if it doesn't yet exist
ptab = table();
if exist(pathsFile,"file")
    ptab = readtable(pathsFile);
end

if size(ptab,1)==0
    ptab.program="";
    ptab.path="";
end

for i=1:length(programs)
    if ~ismember(programs(i), ptab.program)
        this_program=programs(i);
        this_path=string(uigetdir(sctoolboxpath, this_program + " path"));
        ptab{end+1,:} = [this_program, this_path];
    end
end

ptab(ptab.program=="",:)=[];

writetable(ptab,pathsFile)