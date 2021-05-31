function [result, Mdl] = regressOut(cells, Y, params)
% returns the residuals of general linear regression of responses Y (count
% matrix) on predictors X. 

%TODO: more params (eg. clipval, which residuals...)
%TODO: bayesian version, with shrinkage toward global?

X = cells(:,params.vars);

clipval=sqrt(size(Y,1));

result = zeros(size(Y));
for i=1:size(Y,2) %loop over genes
    X.y = Y(:,i);
    Mdl = fitlm(X);
    result(:,i) = Mdl.Residuals.Raw;
%     result(:,i) = Mdl.Residuals.Pearson;
%     result(:,i) = Mdl.Residuals.Standardized;
%     result(:,i) = Mdl.Residuals.Studentized;
    
%     Res = table2array(Mdl.Residuals);
%     boxplot(Res)
%     drawnow

    %progress
    if mod(i,floor(size(Y,2)/10))==0
        fprintf('.')
    end
end
% any(result(:)>clipval)
result(result>clipval)=clipval;
result(result<-clipval)=-clipval;

fprintf('\n')


%     Mdl = fitlm(X, Y(:,i),'RobustOpts','on');

% beta = mvregress(X, Y);
% result = Y - [ones(size(X,1),1), X]*beta;