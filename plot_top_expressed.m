function plot_top_expressed(counts,genes,N, ax, excluded_genes)
msz=1;

%another option:
% N = find(cum_mpg>total_frac,1,'first'); %total_frac is parameter, eg .25

n_UMI=sum(counts,2);

if exist('excluded_genes','var') && ~isempty(excluded_genes)
%     excluded_genes=startsWith(genes.name,'Mt-')|startsWith(genes.name,'Rpl')|startsWith(genes.name,'Rps');
    counts(excluded_genes,:)=[];
    n_UMI(excluded_genes,:)=[];
    genes(excluded_genes,:)=[];
end

totalUMI=sum(n_UMI);

if isnumeric(N) && isscalar(N)
    [mpgs,ixs]=sort(n_UMI,'descend');
    gix=ixs(1:N);
    C=full(counts(gix,:));
    cum_mpg=cumsum(mpgs/totalUMI);
    Nprct=cum_mpg(N)*100;
    
else %interpret N as gene list???
    gix=getGeneIndices(N,genes.name);
    mpg=n_UMI(gix);
    Nprct=sum(mpg)/totalUMI*100;
    C=full(counts(gix,:));
    N=length(gix);
end
title_text=num2str(N) + " genes: " + num2str(Nprct,'%.3g') + "% of total counts";

ONES=ones(1,size(C,2));

if ~exist('ax','var')||isempty(ax)
    ax=gca;
end

hold on
for i=1:N
    y=i*ONES+0.6*rand(1,size(C,2))-0.3;
    scatter(C(i,:),y,msz,'k.')
end
ax.YDir='reverse';
ax.YTick=1:N;
ax.YTickLabel=genes.name(gix);
axis tight
xlabel('counts per cell')
title(title_text)