function [selx,sely,indexnr]=lasso(coords, ax)
% modified to just apply to an existing axis. assumes stuff is plotted
% there

% function [selx,sely,indexnr]=lasso(x,y)
% lasso -  enables the selection/encircling of (clusters of) events in a scatter plot by hand 
%          using the mouse
% 
% Input:    x,y                 - a set of points in 2 column vectors.
% Output:   selx,sely,indexnr   - a set of selected points in 3 column vectors 
% 
% Note:   After the scatter plot is given, selection by mouse is started after any key press. 
%         This is done to be able to ZOOM or CHANGE AXES etc. in the representation before selection 
%         by mouse.
%         Encircling is done by pressing subsequently the LEFT button mouse at the requested positions 
%         in a scatter plot.
%         Closing the loop is done by a RIGHT button press.
%         
% T.Rutten V2.0/9/2003
% plot(x,y,'.')

x=coords(:,1);
y=coords(:,2);

if ~exist('ax','var')
    ax=gca;
end
axes(ax)

las_x=[];
las_y=[];
c=1;
key=0;
disp('press a KEY to start selection by mouse, LEFT mouse button for selection, RIGHT button closes loop')
while key==0
    key=waitforbuttonpress;
%     pause(0.2) %why?
end

hl1=[];
while c==1 
    [a,b,c]=ginput(1);
    delete(hl1)
    las_x=[las_x;a];
    las_y=[las_y;b];
    hl1=line(las_x,las_y);
end


las_x(length(las_x)+1)=las_x(1);
las_y(length(las_y)+1)=las_y(1);
delete(hl1)
line(las_x,las_y,ones(size(las_x)))  %z-coordinate to put it on top

% pause(.2) %why?

in=inpolygon(x,y,las_x,las_y);
indexnr=find(in>0);
selx=x(indexnr);
sely=y(indexnr);

% figure,
line(selx,sely,ones(size(selx)),'color','g','marker','.','linestyle','none');

% legend(num2str([length(x)-length(selx);length(selx)]));
