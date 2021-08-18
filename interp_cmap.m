function cmap=interp_cmap(c1,cmid,c2,n1,n2)
arguments
    c1=[0.05,0.19,0.42]
    cmid=[0.75,0.75,0.75]
    c2=[0.41,0,0.05]
    n1=128
    n2=128
end
%interpolated colormap between c1-cmid-c2
R=[linspace(c1(1),cmid(1),n1),cmid(1),linspace(cmid(1),c2(1),n2)];
G=[linspace(c1(2),cmid(2),n1),cmid(2),linspace(cmid(2),c2(2),n2)];
B=[linspace(c1(3),cmid(3),n1),cmid(3),linspace(cmid(3),c2(3),n2)];

cmap=[R(:),G(:),B(:)];
