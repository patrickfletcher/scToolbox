function plotPairwiseHeatmap(TF, names, fh)
X=TF*TF';
X=X-diag(diag(X));
X(X==0)=nan;

if exist('fh','var')
    figure(fh)
else
    fh=figure();clf
end

hhm=heatmap(X);
hhm.YDisplayLabels=names;
hhm.XDisplayLabels=names;
hhm.MissingDataColor='w';
% title('Ambiguous')
colorbar off
ax=gca;
ax.ColorLimits=[0,max(X(:))/3];
