function [thr,B]=geneThresholdOtsu(X,nz)
%returns threshold for histogram using Otsu's method, one per row of X

levels=256;

if ~exist('nz','var')
    nz=false;
end

% if ~any(X(:)<0)
% ix=find(sum(X,2)>0)'; %work only on rows with at least one non-zero

nGenes=length(X(:,1));
B=zeros(nGenes,levels);
thr=zeros(nGenes,1);
for i=1:size(X,1)
    edges=linspace(min(X(i,:)),max(X(i,:)),levels+1);
    thisX=X(i,:);
    if nz
        thisX=thisX(thisX>0);
    end
    [histcts,edges]=histcounts(thisX,edges);
    [thrBin,btw]=otsu(histcts,levels);
    nextLowerNonEmptyBin=find(histcts(1:thrBin)>0,1,'last');
    thr(i,1)=sum(edges([nextLowerNonEmptyBin+1,thrBin+1]))/2; %average of lower edge of upper bin and upper edge of lower bin
    B(i,:)=btw;
end

end

%TODO: vectorize?
%TODO: write in terminology of Otsu's original paper (1979)
function [level,btw] = otsu(histogramCounts,levels)
    level=1;
    total = sum(histogramCounts); % '''total''' is the number of pixels in the given image. 
    % OTSU automatic thresholding method
    sumB = 0;
    wB = 0;
    maximum = 0.0;
    sum1 = dot( (0:levels-1), histogramCounts); 
    btw=zeros(1,levels);
    for ii=1:levels
        wB = wB + histogramCounts(ii);
        wF = total - wB;
        if (wB == 0 || wF == 0)
            continue;
        end
        sumB = sumB +  (ii-1) * histogramCounts(ii);
        mF = (sum1 - sumB) / wF;
        between = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
        if ( between >= maximum )
            level = ii;
            maximum = between;
        end
        btw(ii)=between;
    end
end



% % OTSU automatic thresholding method
%     level=1;
%     total = sum(histogramCounts); %the number of pixels in the given image. (1 in Otsu79)
%     muK = 0;
%     wK = 0;
%     maximum = 0.0;
%     muT = dot( (0:levels-1), histogramCounts); 
%     between=zeros(1,levels);
%     for i=1:levels
%         wK = wK + histogramCounts(i);
%         w1 = total - wK;
%         if (wK == 0 || w1 == 0)
%             continue;
%         end
%         muK = muK +  (i-1) * histogramCounts(i); 
%         m0 = (muT - muK) / w1; %w0*mu0=w(k)*(mu(k)/w(k))=mu(k)
%         sigmaB2 = wK * w1 * ((muK / wK) - m0) * ((muK / wK) - m0);
%         if ( sigmaB2 >= maximum )
%             level = i;
%             maximum = sigmaB2;
%         end
%         between(i)=sigmaB2;
%     end
% end
