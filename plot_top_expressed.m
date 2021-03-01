function plot_top_expressed(counts,genes,N, figID, no_house_keeping)
msz=1;

%another option:
% N = find(cum_mpg>total_frac,1,'first'); %total_frac is parameter, eg .25

if exist('no_house_keeping','var') && no_house_keeping
    rm_genes=startsWith(genes.name,'Mt-')|startsWith(genes.name,'Rpl')|startsWith(genes.name,'Rps');
    counts(rm_genes,:)=[];
    genes(rm_genes,:)=[];
end
    
if isnumeric(N) && isscalar(N)
    [mpgs,ixs]=sort(genes.n_umi,'descend');
    gix=ixs(1:N);
    C=counts(gix,:);
    cum_mpg=cumsum(mpgs/sum(genes.n_umi));
    Nprct=cum_mpg(N)*100;
    
else %interpret N as gene list???
    gix=getGeneIndices(N,genes.name);
    mpg=genes.n_umi(gix);
    Nprct=sum(mpg)/sum(genes.n_umi)*100;
    C=counts(gix,:);
    N=length(gix);
end
title_text=num2str(N) + " genes: " + num2str(Nprct,'%.3g') + "% of total counts";

ONES=ones(1,size(C,2));

figure(figID);clf
ax=axes();
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