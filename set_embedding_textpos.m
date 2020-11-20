function [textpos,textlab] = set_embedding_textpos(COORDS, ident, cols, dopause)

if ~exist('dopause','var'), dopause=true; end

catnames=categories(ident);
ctcounts=countcats(ident);
catnames(ctcounts==0)=[];
catnames(catnames=="Unc")=[]; %usually not a distinct cluster
catnames(catnames=="Amb")=[]; %usually not a distinct cluster
textpos=zeros(length(catnames),2);
for i=1:length(catnames)
    textpos(i,:)=mean(COORDS(ident==catnames{i},:));
end
textpos=[textpos,ones(length(catnames),2)];

% if ~exist('textlab','var') 
textlab=catnames;
% end

fh=gcf;

ax=plotScatter(COORDS, 'group', ident, cols, fh);

ht1=text(ax,textpos(:,1),textpos(:,2),ones(size(textpos,1),1),textlab,...
    'HorizontalAlignment','center');
drawnow

%pause to adjust label positions
if dopause
    pause
end

%extract the new positions after the pause
textpos=cat(1,ht1(:).Position);