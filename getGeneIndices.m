function [index, queryGenes, dupeix, notfoundix]=getGeneIndices(queryGenes,geneList,options)
arguments
    queryGenes
    geneList
    options.removeDupes=false
    options.verbose=true
end
%attempt to get querygenes indices from genelist.
% If querygene not present, removes.
% If present more than once, add repeats to foundGenes

%ismember should work too: [lia,lib]=ismember(queryGenes,geneList) -> lia is mask on queryGenes, lib is index into
%geneList. However, lib only gives lowest index

if ischar(queryGenes) %single query gene
    queryGenes={queryGenes};
end
nQuery=length(queryGenes);

dupeix=[];
notfoundix=[];

[foundGenes,ixa,ixb]=intersect(queryGenes,geneList,'stable'); %intersect removes duplicates

index=zeros(nQuery,1);
dupeix=[];
if length(ixb)~=nQuery
    %could have duplicates or query is not in geneList.
    notpresent=setxor(1:length(queryGenes),ixa);
    isdupe=ismember(queryGenes(notpresent),foundGenes);
    if any(isdupe)
        dupeix=notpresent(isdupe);
%         dupes=strcat(queryGenes(notpresent(isdupe)), " ");
%         disp(['duplicate query genes: ', dupes{:} ]);
    end
    if any(~isdupe)
        notfoundix=notpresent(~isdupe);
        notfound=strcat(queryGenes(notpresent(~isdupe)), " ");
        if options.verbose
            disp(['not found: ', notfound{:} ]);
        end
    end
end

index(setxor(1:nQuery,[dupeix;notfoundix],'stable'))=ixb;

if options.removeDupes
    index(dupeix)=[];
    queryGenes(dupeix)=[];
else
    for i=1:length(dupeix)
        index(dupeix(i))=index(find(strcmp(queryGenes,queryGenes(dupeix(i))),1,'first'));
    end
end

index(notfoundix)=[];
queryGenes(notfoundix)=[];

%old way, very slow!
% if ~iscell(queryGenes)
%     queryGenes={queryGenes};
% end
% index=[];
% foundGenes=[];
% 
% %do one gene at a time to maintain the query order
% for i=1:length(queryGenes)
%     thisQuery=queryGenes{i};
%     if ~iscell(thisQuery)
%         thisQuery={thisQuery};
%     end
% %     if isscalar(thisQuery)
% %         thisIx=find(strcmp(geneList,thisQuery));
% %         index=[index;thisIx];
% %         foundGenes=[foundGenes;geneList(thisIx)];
% %     else
%         for j=1:length(thisQuery)
%             thisIx=find(strcmp(geneList,thisQuery{j}));
%             index=[index;thisIx];
%             foundGenes=[foundGenes;geneList(thisIx)];
%             
%             if isempty(thisIx)
%                 disp([thisQuery{j} ' was not found']);
%             elseif numel(thisIx)>1
%                 disp([thisQuery{j} ' has ' num2str(numel(thisIx)) ' separate gene IDs']);
%             end
%             
%         end
% %     end
% end