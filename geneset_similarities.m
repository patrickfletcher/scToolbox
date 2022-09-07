function Sim = geneset_similarities(items, groups, metric,options)
arguments
    items
    groups
    metric = "jaccard"
    options.do_plot
end

switch lower(metric)
    case 'jaccard'
        simfun=@(A,B) nnz(A&B)/(nnz(A|B));
    case 'overlap'
        simfun=@(A,B) nnz(A&B)/min(nnz(A),nnz(B));
    case 'dice'
        simfun=@(A,B) 2*nnz(A&B)/(nnz(A)+nnz(B));
    case 'correlation'
        simfun=@(A,B) corr(A(:),B(:));
    case 'nmi'
        simfun=@(A,B) nmi(A,B);
end
       
%one-hot encode gene lists:
uitems=unique(items);
ugroups=unique(groups,'stable');
K=length(ugroups);
I=false(length(uitems),K);
for i=1:K
    thisitems=items(groups==ugroups{i});
    I(ismember(uitems,thisitems),i)=true;
end
 
% now get the similarity matrix
Sim = zeros(K);
for i=1:K
    A=I(:,i);
    for j=1:K
        B=I(:,j);
        Sim(i,j)=simfun(A,B);
    end
end

% Sim = zeros(K);
% for i=1:K-1
%     A=I(i,:);
%     for j=i+1:K
%         B=I(j,:);
%         Sim(i,j)=simfun(A,B);
%     end
% end

if options.do_plot
D=pdist(Sim);
Z=linkage(D,'complete');
o=optimalleaforder(Z,D);

figure()
imagesc(Sim(o,o))
xticks(1:length(ugroups))
xticklabels(ugroups(o))
yticks(1:length(ugroups))
yticklabels(ugroups(o))
axis square
end
