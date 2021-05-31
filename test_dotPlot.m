%plot tests

fh1=figure(1);clf
fh1.Units='inches';
fh1.Position=[2,1,4,2];

sp_params.marg_h=[0.2,0.35];
sp_params.marg_w=[0.1,0.05];
sp_params.minArea=1;
sp_params.maxArea=12;
sp_params.prct_leg=[5,50,100]; %# of entries in legend
sp_params.min_prct=5;
sp_params.cb_prct_gap=0.1;
sp_params.cb_gap=0.025;
sp_params.cb_width=0.05;
sp_params.prct_leg_height=0.25;
sp_params.cblabel='log_{10} mean expression';

% PRCT=(0:5:100)'*(0:5:100)/100;
PRCT=0:5:100;
ctnames=string(PRCT);
EXPR=[linspace(0,3,size(PRCT,2));ones(size(PRCT))];
PRCT=[50*ones(size(PRCT));PRCT];

Glab=strcat('g',string(1:size(PRCT,1)))';
% ctnames=strcat('c',string(1:size(PRCT,2)));

[ax,hs,cb]=plotDotPlot(Glab,ctnames,PRCT,EXPR,fh1,[],'none',sp_params);
ax(1).XTickLabelRotation=90;