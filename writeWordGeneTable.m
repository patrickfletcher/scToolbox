function writeWordGeneTable(ActXWord,domGenes,specGenes, groupNames, colors)
% write dominant and specific genes to a word paragraph in open document
% represented by ActXWord
%
% domGenes/specGenes: cell array of nGroups elements, color of each group
% is set by groupColors.  


nGroups=length(groupNames);
for n=1:nGroups
    
        
    thisDom=sort(domGenes{n});
    thisSpec=sort(specGenes{n});
    ActXWord.Selection.Font.ColorIndex=colors{n};
    
    ActXWord.Selection.ItalicRun;
    isBold=ismember(thisDom,thisSpec);
    for g=1:length(thisDom)
        if isBold(g)
            ActXWord.Selection.BoldRun
            ActXWord.Selection.TypeText(thisDom{g})
            ActXWord.Selection.BoldRun
        else
            ActXWord.Selection.TypeText(thisDom{g})
        end
        if g<length(thisDom)
            ActXWord.Selection.TypeText(', ')
        else
            if n<nGroups
                ActXWord.Selection.TypeText('; ')
            else
                ActXWord.Selection.TypeText('. ')
            end
        end
    end
    ActXWord.Selection.ItalicRun;
    
end