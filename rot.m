function [X, R] = rot(X,theta)
switch size(X,2)
    case 2
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    case 3
        t=theta(1); u=theta(2); v=theta(3);
        Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
        Ry = [cos(u) 0 sin(u); 0 1 0; -sin(u) 0 cos(u)];
        Rz = [cos(v) -sin(v) 0; sin(v) cos(v) 0; 0 0 1];
        R=Rx*Ry*Rz;
end

X=X*R;
% X=(R*X')';

%is it orthonormal?
% norm(R'*R - eye(size(R)))
