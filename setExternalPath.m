function setExternalPath(programs, path, options)
arguments
    programs = ["python";"R"]
    path = []
    options.overwrite = true
end
% create/modify a .env file to store paths to python and R (others?)
% it will be a CSV, for easy use as a table


sctoolboxpath=fileparts(mfilename('fullpath'));
pathsFile = fullfile(sctoolboxpath,'extPaths.csv');

%make the file if it doesn't yet exist
ptab = table();
if exist(pathsFile,"file")
    ptab = readtable(pathsFile,TextType="string");
end

if size(ptab,1)==0
    ptab.program="";
    ptab.path="";
end

% if path not specified, uigetdir. if exist, overwrite?
for i=1:length(programs)
    this_program=programs(i);
    programIsPresent = ismember(this_program, ptab.program);
    if isempty(path)
        this_path=string(uigetdir(sctoolboxpath, this_program + " path"));
    else
        this_path = path(i);
    end
    if programIsPresent
        if options.overwrite
            ptab.path(ptab.program==this_program) = this_path;
        end
    else
        ptab{end+1,:} = [this_program, this_path];
    end
end

ptab(ptab.program=="",:)=[];

ptab

writetable(ptab,pathsFile)