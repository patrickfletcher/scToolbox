function [n_min, min_vals] = find_multimodal(X, minN, genes)

if ~exist('minN','var') || isempty(minN), minN=10; end

N = size(X,1);
n_min=zeros(N,1);
min_vals=cell(N,1);

for i=1:N
    x=X(i,:);
    x(x==0)=[];
    if ~isempty(x) && length(x)>minN && max(x)>min(x)
        
        kfit = fitdist(x(:),'Kernel');
%         kfit = fitdist(x(:),'Kernel','Width',kfit.BandWidth*0.5);
%         lnfit = fitdist(x(:),'Lognormal',
        xx = linspace(min(x),max(x),200);
        yy = kfit.pdf(xx);
        min_x = xx(islocalmin(yy));
        min_y = yy(islocalmin(yy));
        max_y = yy(islocalmax(yy));
%         min_cond = (min_y)./max(max_y) > 0.1;
%         min_cond = max(diff(max_y))>
        min_y = min_y(min_cond);
        min_vals{i} = min_x(min_cond);
        n_min(i)=length(min_vals{i});
        doPlot=true;
        if doPlot && n_min(i)>0
            histogram(x,ceil(sqrt(numel(x))),'Normalization','pdf')
            ax=gca;
            line(xx, yy,'color','k')
            line(min_vals{i}, min_y,'color','r','linestyle','none','marker','s')
            line(genes.thr(i)*[1;1], ax.YLim'*ones(1,2),'color','k')
            genes.name(i)
            genes.n_umi(i)
        end
    end
    if mod(i,round(N/10))==0
        fprintf('.')
    end
end