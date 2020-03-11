function result = doUMAP(score, params, figID, valuenames, colors)

disp('Computing UMAP...')
result = params;

rng(params.rngSeed)

if isempty(params.initY)
    initY=[];

elseif isnumeric(params.initY)&&numel(params.initY)==2
    %initY=[pc1,pc2] indicates index of PCs to use as initial points

    pcix=params.initY;
    signix=sign(pcix);
    pcix=abs(pcix);
    initY=score(:,pcix);
    initY(:,1)=signix(1)*initY(:,1);
    initY(:,2)=signix(2)*initY(:,2);
    
%         initY=score(:,params.initY);
%         initY=initY+ 5*(rand(size(initY))-0.5); %random jitter?

elseif size(params.initY,1)==size(score,1)
    initY=params.initY;

end

%non-standard for run_umap api: 'init', 'spread'
% 'init',initY, 'spread',params.spread, ...    

res=run_umap(score,'n_neighbors',params.n_neighbors,...
    'n_epochs',params.n_epochs,'verbose','text','min_dist',params.min_dist,...
    'randomize',params.randomize,'method',params.method,...
    'python',params.python);

result.coords=res;

if exist('figID','var')
    figure(figID);clf
    plotScatter(res,'value',valuenames,colors,figID);
    axis tight
    drawnow
end
