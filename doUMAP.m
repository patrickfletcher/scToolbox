function result = doUMAP(score, params, figID, valuenames, colors)

disp('Computing UMAP...')
result = params;

Pypath = py.sys.path;
MLpath=string(path).split(';');
sctoolpath=MLpath(contains(MLpath,'scToolbox'));
if count(Pypath,sctoolpath) == 0 && ~isempty(sctoolpath)
    insert(Pypath,int32(0),sctoolpath);
end

umap=py.importlib.import_module('umap');

rng(params.rngSeed)

initY=params.initY;
if isempty(initY)
    initY="spectral";

elseif isnumeric(initY)&&numel(initY)==2
    %initY=[pc1,pc2] indicates index of PCs to use as initial points

    pcix=initY;
    signix=sign(pcix);
    pcix=abs(pcix);
    initY=score(:,pcix);
    initY(:,1)=signix(1)*initY(:,1);
    initY(:,2)=signix(2)*initY(:,2);
    initY=py.numpy.array(initY);

elseif size(initY,1)==size(score,1)
    initY=py.numpy.array(initY);

end

%non-standard for run_umap api: 'init', 'spread'
% 'init',initY, 'spread',params.spread, ...    

% res=run_umap(score,'n_neighbors',params.n_neighbors,...
%     'n_epochs',params.n_epochs,'verbose','text','min_dist',params.min_dist,...
%     'randomize',params.randomize,'method',params.method,...
%     'python',params.python);

% kwargs =py.dict(pyargs('n_neighbors',int32(params.n_neighbors),'min_dist',params.min_dist,...
%                 'n_epochs',int32(params.n_epochs),'verbose','text',...
%                 'n_components',int32(2),'verbose',true,'init',initY));
% res=py.runumap.run(score,kwargs);

kwargs = pyargs('n_neighbors',int32(params.n_neighbors),'min_dist',params.min_dist,...
                'n_epochs',int32(params.n_epochs),'verbose','text',...
                'n_components',int32(2),'verbose',true,'init',initY);
res=umap.UMAP(kwargs).fit_transform(score);
res=double(res);

result.coords=res;

if exist('figID','var')
    figure(figID);clf
    plotScatter(res,'value',valuenames,colors,figID);
    axis tight
    drawnow
end
