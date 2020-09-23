function [P, flag] = diptest(X, minN, nboot)

if ~exist('minN','var') || isempty(minN), minN=4; end
if ~exist('nboot','var') || isempty(nboot), nboot=100; end

P=ones(size(X,1),1);
flag(size(X,1),1)="";
% hasNonZeros=sum(X,2)>0;
for i=1:length(P)
    x=X(i,:);
    x(x==0)=[];
    if ~isempty(x) && length(x)>minN && max(x)>min(x)
        [~,P(i)] = hartigansdipsigniftest(x,nboot);
    else
        % set flag to indicate faults
        if nargout==2
            if isempty(x)
                flag(i)="zeros";
            elseif length(x)<=minN
                flag(i)="minN";
            elseif max(x)==min(x)
                flag(i)="identical";
            end
        end
    end
    if mod(i,round(length(P)/10))==0
        fprintf('.')
    end
end