function result = regressOut(X, Y)
% returns the residuals of general linear regression of responses Y (count
% matrix) on predictors X. 

% beta = mvregress(X, Y);
% result = Y - [ones(size(X,1),1), X]*beta;

result = zeros(size(Y));
for i=1:size(Y,2) %loop over genes
    Mdl = fitlm(X, Y(:,i));
    result(:,i) = Mdl.Residuals.Raw;
%     result(:,i) = Mdl.Residuals.Pearson;
%     result(:,i) = Mdl.Residuals.Standardized;
    
    %progress
    if mod(i,floor(size(Y,2))/10)==0
        i
    end
end