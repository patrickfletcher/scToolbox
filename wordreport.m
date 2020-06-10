function wr = wordreport(varargin)
% wordreport: generate a Microsoft Office Word report
% usage: wr = wordreport;           % use default options
% usage: wr = wordreport(filename); % provide filename for report
%
% WORDREPORT creates or opens an existing Microsoft Office Word report and
% provides helper functions to add some content: text, figures, Simulink
% models, Stateflow charts and much more. It also helps in adding or
% updating the table of contents, in setting page numbering or in finding
% text. Actually, it is possible to mimic almost everything you can do
% manually. Just record a macro in Word and analyze the generated VBA code
% to find out how to use the ActiveX technology. Inside MATLAB, you can
% get a list of available properties with instructions like:
%   hdlActiveX = wr.getactivexhdl();
%   get(hdlActiveX);
%   get(hdlActiveX.Selection);
%   invoke(hdlActiveX.Selection);
%
% Returned is a structure containing function handles that enable to
% manipulate the Word document in an object-oriented way.
%
%
% Arguments: (input)
%  filename - (OPTIONAL) - string - Filename of the document to create or
%        open. The user is responsible for checking file existence before
%        trying to create or open it. If file does not exist, it is
%        created.
%
%        DEFAULT: a unique file name generated from current date and time.
%
% Arguments: (output)
%  wr - structure of function handles - this structure is used to mimic an
%        object-oriented programming syntax. For example, the following
%        syntax can be used to add some text: wr.addtext('Some text')
%
%
% Example:
%  Create a new document called 'Foo.doc' and add some content (headings,
%  figures, page breaks, page numbers, table of contents)
%
%     reportFilename = fullfile(pwd,'foo.doc');
%     wr = wordreport(reportFilename);
%     % ---
%     wr.setstyle('Heading 1');
%     wr.addtext('TOC', [1 1]); % line break before and after text
%     wr.createtoc(1, 3);
%     wr.addpagebreak();
%     % ---
%     wr.setstyle('Heading 1');
%     wr.addtext('MATLAB data', [1 1]); % line break before and after text
%     % ---
%     wr.setstyle('Heading 2');
%     wr.addtext('Sample table', [0 1]); % line break after text
%     dataCell = { ...
%         'Test 1', num2str(0.3) , 'OK'; ...
%         'Test 2', num2str(1.8) , 'KO'};
%     [nbRows, nbCols] = size(dataCell);
%     wr.addtable(nbRows, nbCols, dataCell, [1 1]); % line break before table
%     % ---
%     wr.setstyle('Heading 2');
%     wr.addtext('Sample figure', [0 1]); % line break after text
%     figure; plot(1:10);
%     title('Figure 1'); xlabel('Temps [s]'); ylabel('Amplitude [A]');
%     wr.setstyle('Normal');
%     wr.addfigure();
%     % ---
%     wr.addpagenumbers('wdAlignPageNumberRight');
%     wr.updatetoc();
%     % ---
%     wr.close();
%     % ---
%     open(reportFilename);
%
% More examples:
%  It is also possible to add a Simulink model view with ADDMODEL, a
%  Stateflow chart view with ADDSTATEFLOW, a symbol with ADDSYMBOL, ...
%  Refer to details of helper functions below for further information.
%
%
% Notes:
%  When invoked without a semicolon (i.e wr = wordreport(...)), a
%    list of all available helper functions is displayed.
%  When using SETSTYLE, be careful to use a style name matching your
%    language settings. For example: 'Heading 1' for english users but
%    'Titre 1' for french settings.
%
%
% Acknowledgements:
%  WriteToWordFromMatlab - MATLAB Central File Exchange
%    http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9112&objectType=file
%
%
% Author: Laurent Vaylet
% E-mail: laurent.vaylet@gmail.com
% Release: 1.0
% Release date: 12/10/07
% Display some help if neither input nor output is given
if ~nargin && ~nargout
    help wordreport
    return
end
% Default inputs
inArgs = { ...
    [genvarname(['report' datestr(clock, 'yyyymmddHHMMSS')]) '.doc']}; % filename
% Replace default inputs with specified arguments
inArgs(1:nargin) = varargin;
% Various initializations
docFilename  = inArgs{1};
hdlWordDoc   = [];
hdlActiveX   = [];
currentStyle = 'Normal';
% Create and open a new document (or open an existing one)
CreateDoc();
% Assign output argument (structure containing function handles)
wr = struct( ...
    'addtext',        @AddText, ...
    'addsymbol',      @AddSymbol, ...
    'addtable',       @AddTable, ...
    'addpagenumbers', @AddPageNumbers, ...
    'addfigure',      @AddFigure, ...
    'addmodel',       @AddModel, ...
    'addstateflow',   @AddStateflow, ...
    'setstyle',       @SetStyle, ...
    'goto',           @Goto, ...
    'createtoc',      @CreateTOC, ...
    'updatetoc',      @UpdateTOC, ...
    'addpagebreak',   @AddPageBreak, ...
    'findtext',       @FindText, ...
    'select',         @Select, ...
    'getcolumn',      @GetColumn, ...
    'getactivexhdl',  @GetActiveXHandle, ...
    'printmethods',   @PrintMethods, ...
    'close',          @CloseDoc);
% ---------------------------------------
    function CreateDoc
        % CREATEDOC Create a new Word document using ActiveX and save its handle
        % Start an ActiveX session with Word
        hdlActiveX = actxserver('Word.Application');
        hdlActiveX.Visible = true;
        trace(hdlActiveX.Visible);
        if ~exist(docFilename, 'file');
            % Create new document
            hdlWordDoc = invoke(hdlActiveX.Documents, 'Add');
        else
            % Open existing document
            hdlWordDoc = invoke(hdlActiveX.Documents, 'Open', docFilename);
        end
    end % CreateDoc
% ---------------------------------------
    function AddText(varargin)
        % ADDTEXT Add some text to the document, using current style
        % Default inputs
        inArgs = { ...
            '', ...    % no text
            [0 1], ... % line break after text
            [], ...    % automatic color
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        text       = inArgs{1};
        lineBreaks = inArgs{2};
        color      = inArgs{3};
        % Line breaks before
        for k = 1:lineBreaks(1)
            hdlActiveX.Selection.TypeParagraph;
        end
        % Apply specified style to text and insert it
        hdlActiveX.Selection.Style = currentStyle;
        if ~isempty(color)
            hdlActiveX.Selection.Font.Color = color;
        end
        hdlActiveX.Selection.TypeText(text);
        hdlActiveX.Selection.Font.Color = 'wdColorAutomatic'; % Set back to default color
        % Line breaks after
        for k = 1:lineBreaks(2)
            hdlActiveX.Selection.TypeParagraph;
        end
    end % AddText
% ---------------------------------------
    function AddSymbol(varargin)
        % ADDSYMBOL Add a symbol represented by an integer
        % Integer can be found in the Insert/Symbol menu
        % (176 = degree symbol)
        % Default inputs
        inArgs = { ...
            176, ...    % degree (°)  symbol (random choice)
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        symbol = inArgs{1};
        hdlActiveX.Selection.InsertSymbol(symbol);
    end % AddSymbol
% ---------------------------------------
    function AddTable(varargin)
        % ADDTABLE Add a table to the document
        % Default inputs
        inArgs = { ...
            2, ...                    % rows count
            2, ...                    % columns count
            {'1', '2'; '3', '4'}, ... % data
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        nbRows     = inArgs{1};
        nbCols     = inArgs{2};
        dataCell   = inArgs{3};
        lineBreaks = inArgs{4};
        % Line breaks before
        for k = 1:lineBreaks(1)
            hdlActiveX.Selection.TypeParagraph;
        end
        % Create the table
        % Add = handle Add(handle, handle, int32, int32, Variant(Optional))
        hdlActiveX.ActiveDocument.Tables.Add(hdlActiveX.Selection.Range, nbRows, nbCols, 1, 1);
        % Hard-coded optionals
        % first 1 same as DefaultTableBehavior := wdWord9TableBehavior
        % last  1 same as AutoFitBehavior := wdAutoFitContent
        SetStyle('Normal');
        % Write data into table
        for r = 1:nbRows
            for c = 1:nbCols
                % Write data into current cell
                AddText(dataCell{r, c}, [0, 0]);
                if(r*c == nbRows*nbCols)
                    % Done, leave the table
                    hdlActiveX.Selection.MoveDown;
                else % Move on to next cell
                    hdlActiveX.Selection.MoveRight;
                end
            end
        end
        % Line breaks after
        for k = 1:lineBreaks(2)
            hdlActiveX.Selection.TypeParagraph;
        end
    end
% ---------------------------------------
    function AddPageNumbers(varargin)
        % ADDPAGENUMBERS Add page numbering to the document
        % Default inputs
        inArgs = { ...
            'wdAlignPageNumberRight', ... % page number alignment
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        alignStyle = inArgs{1};
        % Make sure the window isn't split
        if (~strcmp(hdlActiveX.ActiveWindow.View.SplitSpecial, 'wdPaneNone'))
            hdlActiveX.Panes(2).Close;
        end
        % Make sure we are in printview
        if (strcmp(hdlActiveX.ActiveWindow.ActivePane.View.Type, 'wdNormalView') || ...
                strcmp(hdlActiveX.ActiveWindow.ActivePane.View.Type, 'wdOutlineView'))
            hdlActiveX.ActiveWindow.ActivePane.View.Type  = 'wdPrintView';
        end
        % Switch to header-footer view
        hdlActiveX.ActiveWindow.ActivePane.View.SeekView = 'wdSeekCurrentPageHeader';
        if hdlActiveX.Selection.HeaderFooter.IsHeader
            hdlActiveX.ActiveWindow.ActivePane.View.SeekView = 'wdSeekCurrentPageFooter';
        else
            hdlActiveX.ActiveWindow.ActivePane.View.SeekView = 'wdSeekCurrentPageHeader';
        end
        % Add page numbers
        % 0 -> don't display on first page
        hdlActiveX.Selection.HeaderFooter.PageNumbers.Add(alignStyle, 0);
        % Switch back to main document view
        hdlActiveX.ActiveWindow.ActivePane.View.SeekView = 'wdSeekMainDocument';
    end % AddPageNumbers
% ---------------------------------------
    function AddFigure
        % ADDFIGURE Add a figure to the document
        % Capture current figure into clipboard
        print -dmeta
        % Paste clipboard content
        invoke(hdlActiveX.Selection, 'Paste'); % or Paste(hdlActiveX.Selection)
%         hdlActiveX.Selection.TypeParagraph; % line break
    end % AddFigure
% ---------------------------------------
    function AddModel(varargin)
        % ADDMODEL Add a Simulink model to the document
        % Default inputs
        inArgs = { ...
            get_param(gcs, 'Name'), ... % model name
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        modelName = inArgs{1};
        % Capture model into clipboard
        try
            print(['-s' modelName], '-dmeta');
        catch
            lasterr
            warning('WordReport:modelCaptureFailed', 'Cannot capture model (check it is open and given name is correct)');
        end
        % Paste clipboard content
        invoke(hdlActiveX.Selection, 'Paste'); % or Paste(hdlActiveX.Selection)
        hdlActiveX.Selection.TypeParagraph; % line break
    end % AddModel
% ---------------------------------------
    function AddStateflow(varargin)
        % ADDSTATEFLOW Add a Stateflow chart to the document
        error(nargchk(1, 1, nargin));
        stateflowName = varargin{1};
        % Capture model into clipboard
        try
            sfprint(stateflowName, 'meta');
        catch
            warning('WordReport:modelCaptureFailed', 'Cannot capture model (check it is open and given name is correct)');
        end
        % Paste clipboard content
        invoke(hdlActiveX.Selection, 'Paste'); % or Paste(hdlActiveX.Selection)
        hdlActiveX.Selection.TypeParagraph; % line break
    end % AddStateflow
% ---------------------------------------
    function SetStyle(varargin)
        % SETSTYLE Set current text style, used later by AddText
        % Default inputs
        inArgs = { ...
            'Normal', ... % default style
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        style = inArgs{1};
        currentStyle = style;
    end % SetStyle
% ---------------------------------------
    function Goto(varargin)
        % GOTO Jump to specified location in document
        % Default inputs
        inArgs = { ...
            'wdGotoBookmark', ... % what
            'wdGotoAbsolute', ... % which
            1, ...                % count
            '', ...               % name
            false};               % delete ?
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        what   = inArgs{1};
        which  = inArgs{2};
        count  = inArgs{3};
        name   = inArgs{4};
        delete = inArgs{5};
        switch what % enum WdGoToItem
            case 'wdGoToBookmark'
                what = -1;
            case 'wdGoToComment'
                what = 6;
            case 'wdGoToEndnote'
                what = 5;
            case 'wdGoToEquation'
                what = 10;
            case 'wdGoToField'
                what = 7;
            case 'wdGoToFootnote'
                what = 4;
            case 'wdGoToGrammaticalError'
                what = 14;
            case 'wdGoToGraphic'
                what = 8;
            case 'wdGoToHeading'
                what = 11;
            case 'wdGoToLine'
                what = 3;
            case 'wdGoToObject'
                what = 9;
            case 'wdGoToPage'
                what = 1;
            case 'wdGoToPercent'
                what = 12;
            case 'wdGoToProofreadingError'
                what = 15;
            case 'wdGoToSection'
                what = 0;
            case 'wdGoToSpellingError'
                what = 13;
            case 'wdGoToTable'
                what = 2;
        end
        switch which % enum WdGoToDirection
            case 'wdGoToAbsolute'
                which = 1;
            case 'wdGoToFirst'
                which = 1;
            case 'wdGoToLast'
                which = -1;
            case 'wdGoToNext'
                which = 2;
            case 'wdGoToPrevious'
                which = 3;
            case 'wdGoToRelative'
                which = 2;
        end
        hdlActiveX.Selection.GoTo(what, which, count, name);
        if delete
            hdlActiveX.Selection.Delete;
        end
    end % Goto
% ---------------------------------------
    function CreateTOC(varargin)
        % CREATETOC Create the table of contents
        % With ActiveDocument
        %     .TablesOfContents.Add Range := Selection.Range, RightAlignPageNumbers :=  _
        %         True, UseHeadingStyles := True, UpperHeadingLevel := 1, _
        %         LowerHeadingLevel := 3, IncludePageNumbers := True, AddedStyles := "", _
        %         UseHyperlinks := True, HidePageNumbersInWeb := True, UseOutlineLevels :=  _
        %         True
        %     .TablesOfContents(1).TabLeader = wdTabLeaderDots
        %     .TablesOfContents.Format = wdIndexIndent
        % End With
        % Default inputs
        inArgs = { ...
            1, ... % upper heading
            3, ... % lower heading
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        upperHeading = inArgs{1};
        lowerHeading = inArgs{2};
        hdlActiveX.ActiveDocument.TablesOfContents.Add(hdlActiveX.Selection.Range, 1, ...
            upperHeading, lowerHeading);
        hdlActiveX.Selection.TypeParagraph; % Line break after TOC
    end % CreateTOC
% ---------------------------------------
    function UpdateTOC
        % UPDATETOC Update the table of contents
        Goto('wdGoToField', 'wdGoToAbsolute', 1, 'TOC', 1); % last 1 to delete the object
        CreateTOC(1, 3);
    end % UpdateTOC
% ---------------------------------------
    function AddPageBreak
        % ADDPAGEBREAK Add a page break to the document
        hdlActiveX.Selection.InsertBreak;
    end % AddPageBreak
% ---------------------------------------
    function found = FindText(varargin)
        % FINDTEXT Find text in document
        % Default inputs
        inArgs = { ...
            '', ...   % text to find
            true, ... % forward (true) or backward (false)
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        textToFind = inArgs{1};
        forward = inArgs{2};
        hdlActiveX.Selection.Find.ClearFormatting;
        hdlActiveX.Selection.Find.Text = textToFind;
        hdlActiveX.Selection.Find.Replacement.Text = '';
        hdlActiveX.Selection.Find.Forward = forward;
        found = hdlActiveX.Selection.Find.Execute;
    end % FindText
% ---------------------------------------
    function Select(varargin)
        % SELECtT Extends or move selection to the left or the right
        % Default inputs
        inArgs = { ...
            'right', ...       % direction
            'wdCharacter', ... % unit (wdCell, wdCharacter, wdWord, ...)
            1, ...             % count
            'wdMove', ...      % extend or move selection ?
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        direction = inArgs{1};
        unit      = inArgs{2};
        count     = inArgs{3};
        extend    = inArgs{4};
        switch unit % enum WdUnits
            case 'wdCharacter'
                unit = 1;
            case 'wdWord'
                unit = 2;
            case 'wdSentence'
                unit = 3;
            case 'wdParagraph'
                unit = 4;
            case 'wdLine'
                unit = 5;
            case 'wdStory'
                unit = 6;
            case 'wdScreen'
                unit = 7;
            case 'wdSection'
                unit = 8;
            case 'wdColumn'
                unit = 9;
            case 'wdRow'
                unit = 10;
            case 'wdWindow'
                unit = 11;
            case 'wdCell'
                unit = 12;
            case 'wdCharacterFormatting'
                unit = 13;
            case 'wdParagraphFormatting'
                unit = 14;
            case 'wdTable'
                unit = 15;
            case 'wdItem'
                unit = 16;
        end
        switch extend % enum WdMovementType
            case 'wdMove'
                extend = 0;
            case 'wdExtend'
                extend = 1;
        end
        switch direction
            case 'left'
                hdlActiveX.Selection.MoveLeft(unit, count, extend);
            case 'right'
                hdlActiveX.Selection.MoveRight(unit, count, extend);
        end
    end % FindText
% ---------------------------------------
    function data = GetColumn(varargin)
        % GETCOLUMN Get data in specified column of current table
        % Default inputs
        inArgs = { ...
            1, ...    % Index of column to capture
            true, ... % Ignore column header ?
            };
        % Replace default inputs with specified arguments
        inArgs(1:nargin) = varargin;
        % Initializations
        idxCol = inArgs{1};
        ignoreHeader = inArgs{2};
        data = [];
        if hdlActiveX.Selection.Information('wdWithInTable') % check cursor is in table
            nbRows = hdlActiveX.Selection.Tables(1).Item(1).Rows.Count - ignoreHeader; % ignoreHeader = true (1) or false (0)
            % Retrieve data in 2nd column
            data = cell(nbRows, 1); % preallocate
            for row = 1:nbRows
                cellText = hdlActiveX.Selection.Tables(1).Item(1).Cell(row+ignoreHeader,idxCol).Range.Text; % row+1 -> ignore header
                data{row,1} = cellText(1:end-2); % end-2 -> ignore line break
            end
        end
    end % GetColumn
% ---------------------------------------
    function res = GetActiveXHandle
        % GETACTIVEXHANDLE Return current document ActiveX handle
        res = hdlActiveX;
    end % GetActiveXHandle
% ---------------------------------------
    function CloseDoc
        % CLOSEDOC Close an open document
        if ~exist(docFilename, 'file')
            % Save file as new
            invoke(hdlWordDoc, 'SaveAs', docFilename, 1);
        else
            % Save existing file
            invoke(hdlWordDoc, 'Save');
        end
        % Close document window
        invoke(hdlWordDoc, 'Close');
        % Quit Word
        invoke(hdlActiveX, 'Quit');
        % Terminate ActiveX
        delete(hdlActiveX);
        hdlWordDoc = [];
        hdlActiveX = [];
    end % CloseDoc
% ---------------------------------------
    function PrintMethods(varargin)
        % PRINTMETHODS Print all available methods for a Word ActiveX
        category = varargin{1};
        headingString = varargin{2};
        SetStyle([headingString '3']);
        text = strcat(category, '-methods');
        AddText(text, [1,1]);
        SetStyle('Normal');
        text = ['Methods called from Matlab as: hdlActiveX.' category '.MethodName(xxx)'];
        AddText(text, [0 0]);
        text = 'Ignore the first parameter "handle". ';
        AddText(text, [1 3]);
        structMethods = invoke(hdlActiveX.(category));
        cellMethods = struct2cell(structMethods);
        for i = 1:length(cellMethods)
            methodString = cellMethods{i};
            AddText(methodString, [0 1]);
        end
    end % PrintMethods
% ---------------------------------------
end % wordreport
