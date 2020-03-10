function [s,U,S,V] = svdecon(X)
% Input:
% X : m x n matrix
%
% Output:
% X = U*S*V'
%
% Description:
% Does equivalent to svd(X,'econ') but faster
%
% Vipin Vijayan (2014)
%
% https://www.mathworks.com/matlabcentral/fileexchange/47132-fast-svd-and-pca
%
% modified to return just s if desired, without forming any matrices
%X = bsxfun(@minus,X,mean(X,2));
[m,n] = size(X);
if  m <= n
%     C = X*X';
    [U,D] = eig(X'*X);
%     clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    s = sqrt(d);  
    
    if nargout > 2
        U = U(:,ix);  
        V = X'*U;
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
%     C = X'*X; 
    [V,D] = eig(X'*X);
%     clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    s = sqrt(d);
    
    if nargout > 2
        V = V(:,ix);    
        U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
        %s = sqrt(sum(U.^2,1))';
        U = bsxfun(@(x,c)x./c, U, s');
        S = diag(s);
    end
end
