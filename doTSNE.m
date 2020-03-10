function result = doTSNE(score, params, figID, valuenames, colors)

disp('Computing tSNE...')
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

res=tsne(score,'perplexity',params.perplexity,'learnRate',params.learnrate,'exaggeration',params.exagg,...
    'initialY',initY, 'verbose',2,'NumPrint',100,'options',statset('maxIter',params.maxiter));

result.coords=res;
    
if exist('figID','var')
    figure(figID);clf
    plotScatter(res,'value',valuenames,colors,figID);
    axis tight
    drawnow
end