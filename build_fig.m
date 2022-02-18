function fig = build_fig(ax,layoutopts,figopts)
arguments
    ax
    layoutopts.grid=[1,1]      %overall grid dimensions
    layoutopts.margins = 0.05  %[left,right,bottom,top]. scalar - all same, 2 vals - [xmarg, ymarg]
    layoutopts.gap=0.05
    figopts.?matlab.ui.Figure  %use these to set all fig properties
end
% BUILD_FIG construct a layout of panels with panel labels from existing
% axes.

%Could be a simple class: set/get functions to apply properties to all
%panels?

% use cases: 
% - list of axes, global panel grid dimensions, list of panel spans, global margins plus gaps
% - array of struct

%two types of info needed:
% 1 - per panel info: ax, tilespan, panel annotations, etc
% 2 - global figure info: figure properties, overall title, global margins/gaps 

%panel info:
% - ax to copy
% - label
% - index, span into grid

%how to handle legends?

%figopts:
% - need [w,h] + units, maybe [x,y,w,h]?  figpos=[x,y,w,h]
%set some defaults for figopts?
if ~isfield(figopts,"Units")
end


%make a new figure

%compute a grid for panels and add margins/gaps (from tight_supblot)
nr=2;
nc=2;
marg_w=[0.05,0.05];
marg_h=[0.05,0.05];
gap=[0.05,0.05];

axw = (1-sum(marg_w)-(nc-1)*gap(2))/nc;
axh = (1-sum(marg_h)-(nr-1)*gap(1))/nr; 
py = 1-marg_h(2)-axh; 
pos=zeros(nr*nc,4);
for i = 1:nr
    px = marg_w(1);
    for j = 1:nc
        if py<0, py=0; end
        ix=nc*(i-1)+j;
        pos(ix,:)=[px py axw axh];
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end

%copyobj to duplicate input axes into their new axes in the layout

%add annotations: panel labels, titles, etc.