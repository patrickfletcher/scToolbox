function plot_excl_inter(excl_inter, labels)

counts=cellfun(@length,excl_inter);
zerocounts=counts==0;

if length(labels)<=3
   
% ax=tight_subplot(1,1,1,0.05,0.05,0.05);
counts_eps=counts;
counts_eps(zerocounts)=0.02*sum(counts);
[h,S]=venn(counts_eps,'FaceAlpha', 0.6);
axis equal

axis off
% for i=1:length(h)
%     h(i).FaceColor=cols(i,:);
% end

for i = 1:length(counts)
  text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(counts(i)));
end

% legend(labels)
% legend boxoff
% drawnow

else
% if more than 3, use an upset plot

upsetplot(setsD,labels)
end
