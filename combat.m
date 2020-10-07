function result = combat(X, group)

% scToolbox needs to be on Python path:
Pypath = py.sys.path;
MLpath=string(path).split(';');
sctoolpath=MLpath(contains(MLpath,'scToolbox'));
if count(Pypath,sctoolpath) == 0 && ~isempty(sctoolpath)
    insert(Pypath,int32(0),sctoolpath);
end

% import/reload (not sure this is necessary)
combat = py.importlib.import_module('combat');
py.importlib.reload(combat)


Xpy = py.numpy.array(X,'float32');
Gpy = int64(group); %grouping var as int64 (can't seem to pass other)
Xcb = combat.combat(Xpy, Gpy);
% Xcb.head()

result = double(Xcb.to_numpy());