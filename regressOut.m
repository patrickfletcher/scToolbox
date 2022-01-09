function [result, Mdl] = regressOut(cellTab, Y, vars_to_regress, params)
arguments
    cellTab
    Y
    vars_to_regress
    params.residual{mustBeMember(params.residual,["Raw","Pearson","Standardized","Studentized"])}='Raw'
    params.clipval=sqrt(size(Y,1))
end
% returns the residuals of general linear regression of responses Y (count
% matrix) on predictors X. 

%TODO: more params (eg. clipval, which residuals...)
%TODO: bayesian version, with shrinkage toward global?
%TODO: fitlme? for grouping var support
%TODO: reuse computations done uniquely on predictors - QR decomp

nObs=height(cellTab);

[m,n]=size(Y);
if m~=nObs && n==nObs
    Y=Y';
end

%rows are observations
X = cellTab{:,vars_to_regress};

% X = QR
% ||y-Xb||^2 = 

% [Q,R]=qr(X,)
% dX = decomposition(X); 

result = zeros(size(Y));
for i=1:size(Y,2) %loop over genes
    Mdl = fitlm(X, Y(:,i));

    result(:,i) = Mdl.Residuals.(params.residual);
    
%     Res = table2array(Mdl.Residuals);
%     boxplot(Res)
%     drawnow

    %progress
    if mod(i,floor(size(Y,2)/10))==0
        fprintf('.')
    end
end
% any(result(:)>clipval)
result(result>params.clipval)=params.clipval;
result(result<-params.clipval)=-params.clipval;

fprintf('\n')