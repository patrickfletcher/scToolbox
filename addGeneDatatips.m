function txt = addGeneDatatips(~,event_obj,geneName)
% Customizes text of data tips:
%    dcm_obj = datacursormode(fig);
%    set(dcm_obj,'UpdateFcn',{@addGeneDatatips,genes.name})

pos = get(event_obj,'Position');
idx = get(event_obj, 'DataIndex');
txt = {['\it ',geneName{idx}],...
       ['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))]};