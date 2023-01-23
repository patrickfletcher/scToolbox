function [hcb, cblab]=makeCB(ax,rectpos,gap,dims,placement,justify, label)
arguments
    ax=[]
    rectpos=[]
    gap=0
    dims=[0.3,0.02]
    placement {mustBeMember(placement,["east","west","north","south"])}='east'
    justify {mustBeMember(justify,["low","midlo","mid","midhi","high"])}='high'
    label = ""
end

if isempty(ax)
    ax=gca();
end
if isempty(rectpos)
    rectpos=ax.Position;
end

%common features of east/west or north/south
switch placement
    case {'east','west'}
        cbh=dims(1)*rectpos(4);
        cbw=dims(2);
        switch justify
            case 'low'
                cby=rectpos(2);
            case 'midlo'
                cby=rectpos(2)+(rectpos(4)-cbh)*0.25;
            case 'mid'
                cby=rectpos(2)+(rectpos(4)-cbh)*0.5;
            case 'midhi'
                cby=rectpos(2)+(rectpos(4)-cbh)*0.75;
            case 'high'
                cby=rectpos(2)+rectpos(4)-cbh;
        end
    case {'north','south'}
        cbw=dims(1)*rectpos(3);
        cbh=dims(2);
        switch justify
            case 'low'
                cbx=rectpos(1);
            case 'midlo'
                cbx=rectpos(1)+(rectpos(3)-cbw)*0.25;
            case 'mid'
                cbx=rectpos(1)+(rectpos(3)-cbw)*0.5;
            case 'midhi'
                cbx=rectpos(1)+(rectpos(3)-cbw)*0.75;
            case 'high'
                cbx=rectpos(1)+rectpos(3)-cbw;
        end
end
%unique features of each side: gap from axis
switch placement
    case 'east'
        cbx=rectpos(1)+rectpos(3)+gap;
    case 'west'
        cbx=rectpos(1)-gap-cbw;
    case 'north'
        cby=rectpos(2)+rectpos(4)+gap;
    case 'south'
        cby=rectpos(2)-gap-cbh;
end

%     [cbx,cby,cbw,cbh]
hcb=colorbar(ax,placement);
hcb.Position=[cbx,cby,cbw,cbh];

cblab = [];
if ~isempty(label)
    switch placement
        case {'east','west'}
            cblab=ylabel(hcb,label);
        case {'north','south'}
            cblab=xlabel(hcb,label);
    end
end