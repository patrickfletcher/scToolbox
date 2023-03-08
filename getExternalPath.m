function this_path = getExternalPath(program)
arguments
    program
end

sctoolboxpath=fileparts(mfilename('fullpath'));
pathsFile = fullfile(sctoolboxpath,'extPaths.csv');

if ~exist(pathsFile,"file")
    setExternalPath(program);
end

ptab = readtable(pathsFile);

this_path = ptab.path{ptab.program==program};