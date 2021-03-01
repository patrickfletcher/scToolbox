function writeWordParagraph(ActXWord,words,isBold,isItalic,colorGroup,colors)
% write word list to a word paragraph in open document
% represented by ActXWord
% each word has toggles for bold/italic and has color given by
% colors(colorGroup(i))

doColors=true;
if ~exist('colorGroup','var')
    doColors=false;
end

boldOn=false;
italicOn=false;
for i=1:length(words)
    if isBold(i) && ~boldOn
        ActXWord.Selection.Font.Bold=1;
        boldOn=true;
    elseif ~isBold(i) && boldOn
        ActXWord.Selection.Font.Bold=0;
        boldOn=false;
    end
    
    if isItalic(i)&&~italicOn
        ActXWord.Selection.Font.Italic=1;
        italicOn=true;
    elseif ~isItalic(i)&&italicOn
        ActXWord.Selection.Font.Italic=0;
        italicOn=false;
    end
    
    if doColors
        ActXWord.Selection.Font.ColorIndex=colors{colorGroup(i)};
%     ActXWord.Selection.Font.TextColor.RGB=
    end

    ActXWord.Selection.TypeText(words{i})
    
    %restore default
%     ActXWord.Selection.Font.Bold=0;
%     ActXWord.Selection.Font.Italic=0;
%     ActXWord.Selection.Font.ColorIndex='wdAuto';
    
    if i<length(words)
        ActXWord.Selection.TypeText(', ')
    else
        ActXWord.Selection.TypeText('. ')
    end
end

ActXWord.Selection.Font.Bold=0;
ActXWord.Selection.Font.Italic=0;
ActXWord.Selection.Font.ColorIndex='wdAuto';