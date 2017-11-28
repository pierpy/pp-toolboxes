% Edit cellarray of mixed numerical and string content (same content per column)
%
% function Tdata = table_dialog(Tdata,Theader,Title,Mode)
% function Tdata = table_dialog(Tdata,Theader,Title)
% function Tdata = table_dialog(Tdata,Theader)
% function Tdata = table_dialog(Tdata)
%
%   Tdata   =   numerical array OR cell array
%   Theader =   cellarray (size = 1 x size(Tdata,2)) with strings
%   Title   =   string
%   Mode    =   0 = no adding of columns or rows allowed
%               1 = adding of rows allowed
%               2 = adding of rows and columns allowed
%                
% Written by F. Hatz, Neurology Basel 2013

function Tdata = lab_table_dialog(Tdata,Theader,Title,Mode,ColumnFormat,FigColor)

if ~exist('Tdata','var')
    disp('  no input')
    return
end
if ~exist('Mode','var')
    Mode = 2;
end
if ~exist('Title','var')
    Title = 'Table';
end
if ~isempty(Tdata) & exist('Theader','var') & ~isempty(Theader)
    if size(Theader,1) > 1 & size(Theader,2) == 1
        Theader = Theader';
    end
    if size(Theader,2) ~= size(Tdata,2)
        clearvars Theader
    end
end
if ~exist('Theader','var') & ~isempty(Tdata);
    Theader = cellstr(num2str((1:size(Tdata,2))'))';
end
if isempty(Tdata) & ~isempty(Theader)
    Tdata = cell(1,size(Theader,2));
end
if exist('ColumnFormat','var') 
    if size(ColumnFormat,1) > 1 & size(ColumnFormat,2) == 1
        ColumnFormat = ColumnFormat';
    end
    if size(ColumnFormat,2) ~= size(Theader,2)
        clearvars ColumnFormat
    end
end

figT = figure;
posF = get(figT, 'position');

% find max length of cells (including Theader)
Cwidth = 0;
for i = 1:length(Theader)
    len = length(Theader{i});
    if(len > Cwidth)
        Cwidth = len;
    end
end
for i=1:size(Tdata,2)
    % Iterate over each row
    for j=1:size(Tdata,1);
        if iscell(Tdata)
            len = length(Tdata{j,i});
        else
            len = length(Tdata(j,i));
        end
        % Store in Cwidth only if len is max length
        if(len > Cwidth)
            Cwidth = len;
        end
    end
end

% Some calibration needed as ColumnWidth is in pixels
Cwidth = Cwidth*10;
if Cwidth < 70
    Cwidth = 70;
elseif Cwidth > 120
    Cwidth = 120;
end
numrows = size(Tdata,1)+1;
if numrows > 25
    numrows = 25;
elseif numrows < 4
    numrows = 4;
end
posF(3) = Cwidth*size(Tdata,2) + 70;
if posF(3) > 1500
    posF(3) = 1500;
    posF(1) = 70;
elseif posF(3) > 600
    posF(1) = 150;
end
posF(4)=(numrows+1)*24+80;
if posF(4) > 800
    posF(4) = 800;
    posF(2) = 70;
elseif posF(4) > 600
    posF(2) = 70;
end
if ~exist('FigColor','var')
    FigColor=get(0,'DefaultUicontrolBackgroundcolor');
end
set(figT,'position',posF,'MenuBar','None','Name',Title,'NumberTitle','off','Color',FigColor);
pos = posF;
pos(1) = 5;
pos(2) = 30;
pos(3) = pos(3) -10;
pos(4) = pos(4) -35;
if ~exist('ColumnFormat','var') 
    ColumnFormat = {};
    for i = 1:size(Tdata,2)
        if iscell(Tdata)
            if isnumeric(Tdata{1,i})
                ColumnFormat = [ColumnFormat {'numeric'}]; %#ok<AGROW>
            elseif islogical(Tdata{1,i})
                ColumnFormat = [ColumnFormat {'logical'}]; %#ok<AGROW>
            else
                ColumnFormat = [ColumnFormat {'char'}]; %#ok<AGROW>
            end
        elseif isnumeric(Tdata(1,i))
            ColumnFormat = [ColumnFormat {'numeric'}]; %#ok<AGROW>
        elseif islogical(Tdata(1,i))
            ColumnFormat = [ColumnFormat {'logical'}]; %#ok<AGROW>
        else
            ColumnFormat = [ColumnFormat {[]}]; %#ok<AGROW>
        end
    end
end
if Mode >= 0
    ColumnEdit = true(1,size(Tdata,2));
    t = uitable(figT,'Data',Tdata,'ColumnName',Theader,'ColumnEdit',ColumnEdit, ...
        'Position',pos,'ColumnWidth',{Cwidth},'ColumnFormat',ColumnFormat, ...
        'Units','Normalized');
else
    ColumnEdit = false(1,size(Tdata,2));
    t = uitable(figT,'Data',Tdata,'ColumnName',[],'RowName',[],'ColumnEdit',ColumnEdit, ...
        'Position',pos,'ColumnWidth',{Cwidth},'ColumnFormat',ColumnFormat, ...
        'Units','Normalized','BackgroundColor',FigColor);
end

tmp = get(t);
if tmp.Extent(3)>tmp.Position(3);
    posF(4) = posF(4)+22;
end
if tmp.Extent(4)>tmp.Position(4);
    posF(3) = posF(3)+22; %#ok<NASGU>
end

% create menu
if Mode >= 0
    J = findjobj(t);
    m1 = uimenu(figT,'Label','File');
    uimenu(m1,'Label','Store','Callback',@(~,~)do_store);
    uimenu(m1,'Label','Import','Callback',@(hObj,evd)tabledialog_import(t));
    uimenu(m1,'Label','Export','Callback',@(hObj,evd)tabledialog_export(t));
    if length(Theader) == 3 & strcmp(Theader,{'name','active','reference'})
        uimenu(m1,'Label','Import Montage','Callback',@(hObj,evd)read_montage(t));
    end
    if size(Tdata,1) == size(Tdata,2)
        uimenu(m1,'Label','Generate','Callback',@(hObj,evd)generate_matrix(t));
    end
    uimenu(m1,'Label','Close','Callback','close;','Separator','on');
    m2 = uimenu(figT,'Label','Edit');
    if Mode == 2
        uimenu(m2,'Label','+ column','Callback',@(hObj,evd)tabledialog_column_add(t,J));
        uimenu(m2,'Label','- column','Callback',@(hObj,evd)tabledialog_column_del(t,J));
    end
    if Mode > 0
        uimenu(m2,'Label','+ row','Separator','on','Callback',@(hObj,evd)tabledialog_row_add(t,J));
        uimenu(m2,'Label','- row','Callback',@(hObj,evd)tabledialog_row_del(t,J));
    end
else
    m1 = uimenu(figT,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Export','Callback',@(hObj,evd)tabledialog_export(t));
    uimenu(m1,'Label','Close','Callback','close;');
end

% create buttons
if Mode >= 0
    uicontrol('style', 'pushbutton', 'units', 'pixels','position', [2,2,45,22],...
        'string','store','fontsize',10,'Callback',@(hObj,evd)do_store,'TooltipString','ok');
    uicontrol('style', 'pushbutton', 'units', 'pixels','position', [52,2,45,22],...
        'string','cancel','fontsize',10,'Callback',@(hObj,evd)do_cancel,'TooltipString','cancel');
end
if Mode > 0
    uicontrol('style', 'pushbutton', 'units', 'pixels','position', [102,2,45,22], ...
        'string','- row','fontsize',10,'Callback',@(hObj,evd)tabledialog_row_del(t,J),'TooltipString','delete row');
    uicontrol('style', 'pushbutton', 'units', 'pixels','position', [152,2,45,22], ...
        'string','+ row','fontsize',10,'Callback',@(hObj,evd)tabledialog_row_add(t,J),'TooltipString','add row');
end
uiwait;

    % start support scripts table_dialog
    function do_store
        Tdata = get(t,'Data');
        close(figT);
    end

    function do_cancel
        Tdata = [];
        close(figT);
    end

    function tabledialog_import(Thandle)
        [filename,filepath] = uigetfile('*.xls;*.xlsx','Select xls-file');
        if filename == 0
            return;
        end
        if filename == 0
            table = [];
        elseif ispc
            [~,~,table] = xlsread(fullfile(filepath,filename));
        else
            [~,~,table] = xlsread(fullfile(filepath,filename),1,'','basic');
        end
        if isempty(table)
            return
        end
        answer = questdlg('First line is titles?','Title settings','Yes','No','No');
        if strcmp(answer,'No')
            titles = cellstr(num2str((1:size(table,2))'));
        else
            titles = table(1,:)';
            table = table(2:end,:);
        end
        Ttmp=get(Thandle);
        Ttmp.ColumnWidth=Ttmp.ColumnWidth(1);
        if Ttmp.ColumnWidth{1} < 70
            Ttmp.ColumnWidth{1} = 70;
        end
        Ttmp.ColumnFormat={};
        for nCol = 1:size(table,2)
            if isnumeric(table{1,nCol})
                Ttmp.ColumnFormat = [Ttmp.ColumnFormat {'numeric'}];
            elseif islogical(table{1,nCol})
                Ttmp.ColumnFormat = [Ttmp.ColumnFormat {'logical'}];
            else
                Ttmp.ColumnFormat = [Ttmp.ColumnFormat {'char'}];
            end
        end
        set(Thandle,'Data',table,'ColumnEdit',true(1,size(table,2)),'ColumnName',titles,'ColumnWidth',Ttmp.ColumnWidth,'ColumnFormat',Ttmp.ColumnFormat);
        Ttmp = get(Thandle);
        Ttmp2 = get(gcf,'Position');
        Ttmp2(3) = Ttmp2(3)*Ttmp.Extent(3)/Ttmp.Position(3);
        Ttmp2(4) = Ttmp2(4)*Ttmp.Extent(4)/Ttmp.Position(4);
        if Ttmp2(3) > 1500
            Ttmp2(3) = 1500;
            Ttmp2(4) = Ttmp2(4)+22;
        end
        if Ttmp2(4) < 150
            Ttmp2(4) = 150;
        elseif Ttmp2(4) > 800
            Ttmp2(4) = 800;
            Ttmp2(2) = 70;
            Ttmp2(3) = Ttmp2(3)+22;
        end
        set(gcf,'Position',Ttmp2);
    end

    function tabledialog_export(Thandle)
        [filename,filepath] = uiputfile('.xls','File to store');
        if filename == 0
            return
        end
        Ttmp=get(Thandle);
        if ~iscell(Ttmp.Data)
            xlsout = cat(1,Ttmp.ColumnName',num2cell(Ttmp.Data));
        else
            xlsout = cat(1,Ttmp.ColumnName',Ttmp.Data);
        end
        lab_write_xls(fullfile(filepath,filename),xlsout);
    end

    function tabledialog_column_add(Thandle,Jhandle)
        cols = Jhandle.getComponent(0).getComponent(0).getSelectedColumns+1;
        Ttmp=get(Thandle);
        if isempty(cols)
            cols = size(Ttmp.Data,2);
        end
        datatmp = Ttmp.Data(:,1:cols);
        headertmp = Ttmp.ColumnName(1:cols);
        formattmp = Ttmp.ColumnFormat(1:cols);
        if cols < size(Ttmp.Data,2)
            datatmp(:,cols+2:size(Ttmp.Data,2)+1) = Ttmp.Data(:,cols+1:end);
            headertmp(cols+2:size(Ttmp.Data,2)+1) = Ttmp.ColumnName(cols+1:end);
            formattmp(cols+2:size(Ttmp.Data,2)+1) = Ttmp.ColumnFormat(cols+1:end);
        end
        Ttmp.Data = datatmp;
        Ttmp.ColumnName = headertmp;
        Ttmp.ColumnFormat = formattmp;
        clearvars datatmp headertmp formattmp
        nCol = cols+1;
        if isnumeric(Ttmp.Data)
            Ttmp.Data(1,nCol)=0;
        elseif strcmp(Ttmp.ColumnFormat{cols},'numeric')
            Ttmp.Data(:,nCol)=repmat({[]},size(Ttmp.Data,1),1);
            Ttmp.ColumnFormat{nCol} = 'numeric';
        elseif strcmp(Ttmp.ColumnFormat{cols},'logical')
            Ttmp.Data(:,nCol)=num2cell(false(size(Ttmp.Data,1),1));
            Ttmp.ColumnFormat{nCol} = 'logical';
        else
            Ttmp.Data(:,nCol)=repmat({' '},size(Ttmp.Data,1),1);
            Ttmp.ColumnFormat{nCol} = 'char';
        end
        Ttmp.ColumnName{nCol}='new';
        Ttmp.ColumnWidth=Ttmp.ColumnWidth(1);
        if Ttmp.ColumnWidth{1} < 70
            Ttmp.ColumnWidth{1}=70;
        end
        set(Thandle,'Data',Ttmp.Data,'ColumnName',Ttmp.ColumnName,'ColumnEdit',true(1,size(Ttmp.Data,2)),'ColumnWidth',Ttmp.ColumnWidth,'ColumnFormat',Ttmp.ColumnFormat);
        Ttmp=get(Thandle);
        Ttmp2=get(gcf,'Position');
        Ttmp2(3)=Ttmp2(3)*Ttmp.Extent(3)/Ttmp.Position(3);
        Ttmp2(4)=Ttmp2(4)*Ttmp.Extent(4)/Ttmp.Position(4);
        if Ttmp2(3) > 1500
            Ttmp2(3) = 1500;
            Ttmp2(4) = Ttmp2(4)+22;
        end
        if Ttmp2(4)<150
            Ttmp2(4)=150;
        elseif Ttmp2(4)>800
            Ttmp2(4)=800;
            Ttmp2(2)=70;
            Ttmp2(3)=Ttmp2(3)+22;
        end
        set(gcf,'Position',Ttmp2);
    end

    function tabledialog_column_del(Thandle,Jhandle)
        cols = Jhandle.getComponent(0).getComponent(0).getSelectedColumns+1;
        Ttmp=get(Thandle);
        if isempty(cols)
            cols = size(Ttmp.Data,2);
        end
        if cols > 1
            datatmp = Ttmp.Data(:,1:cols-1);
            headertmp = Ttmp.ColumnName(1:cols-1);
            formattmp = Ttmp.ColumnFormat(1:cols-1);
        end
        if cols < size(Ttmp.Data,2)
            datatmp(:,cols:size(Ttmp.Data,2)-1) = Ttmp.Data(:,cols+1:end);
            headertmp(cols:size(Ttmp.Data,2)-1) = Ttmp.ColumnName(cols+1:end);
            formattmp(cols:size(Ttmp.Data,2)-1) = Ttmp.ColumnFormat(cols+1:end);
        end
        if ~exist('datatmp','var')
            Ttmp.Data = {[]};
            Ttmp.ColumnName = {'New'};
            Ttmp.ColumnFormat = {'numeric'};
        else
            Ttmp.Data = datatmp;
            Ttmp.ColumnName = headertmp;
            Ttmp.ColumnFormat = formattmp;
        end
        clearvars datatmp headertmp formatstmp
        Ttmp.ColumnWidth=Ttmp.ColumnWidth(1);
        if Ttmp.ColumnWidth{1} < 70
            Ttmp.ColumnWidth{1} = 70;
        end
        set(Thandle,'Data',Ttmp.Data,'ColumnFormat',Ttmp.ColumnFormat,'ColumnName',Ttmp.ColumnName,'ColumnEdit',true(1,size(Ttmp.Data,2)),'ColumnWidth',Ttmp.ColumnWidth);
        Ttmp = get(Thandle);
        Ttmp2 = get(gcf,'Position');
        Ttmp2(3) = Ttmp2(3)*Ttmp.Extent(3)/Ttmp.Position(3);
        Ttmp2(4) = Ttmp2(4)*Ttmp.Extent(4)/Ttmp.Position(4);
        if Ttmp2(3) > 1500
            Ttmp2(3) = 1500;
            Ttmp2(4) = Ttmp2(4)+22;
        end
        if Ttmp2(4) < 150;
            Ttmp2(4) = 150;
        elseif Ttmp2(4) > 800
            Ttmp2(4) = 800;
            Ttmp2(2) = 70;
            Ttmp2(3) = Ttmp2(3)+22;
        end
        set(gcf,'Position',Ttmp2);
    end

    function tabledialog_row_add(Thandle,Jhandle)
        rows = Jhandle.getComponent(0).getComponent(0).getSelectedRows+1;
        Ttmp=get(Thandle);
        if isempty(rows)
            nRow = size(Ttmp.Data,1)+1;
        else
            datatmp(1:rows,:) = Ttmp.Data(1:rows,:);
            if size(Ttmp.Data,1) > rows
                datatmp(rows+2:size(Ttmp.Data,1)+1,:) = Ttmp.Data(rows+1:end,:);
            end
            nRow = rows+1;
            Ttmp.Data = datatmp;
            clearvars datatmp;
        end
        for nCol = 1:size(Ttmp.Data,2);
            if isnumeric(Ttmp.Data);
                Ttmp.Data(nRow,nCol) = 0;
            elseif strcmp(Ttmp.ColumnFormat{nCol},'numeric');
                Ttmp.Data{nRow,nCol} = [];
            elseif strcmp(Ttmp.ColumnFormat{nCol},'logical');
                Ttmp.Data{nRow,nCol} = false;
            else
                Ttmp.Data{nRow,nCol} = ' ';
            end;
        end;
        set(Thandle,'Data',Ttmp.Data);
        Ttmp = get(Thandle);
        Ttmp2 = get(gcf,'Position');
        Ttmp2(3) = Ttmp2(3)*Ttmp.Extent(3)/Ttmp.Position(3);
        Ttmp2(4) = Ttmp2(4)*Ttmp.Extent(4)/Ttmp.Position(4);
        if Ttmp2(3) > 1500
            Ttmp2(3) = 1500;
            Ttmp2(4) = Ttmp2(4)+22;
        end
        if Ttmp2(4) < 150
            Ttmp2(4) = 150;
        elseif Ttmp2(4) > 800
            Ttmp2(4) = 800;
            Ttmp2(2) = 70;
            Ttmp2(3) = Ttmp2(3)+22;
        end
        set(gcf,'Position',Ttmp2);
    end

    function tabledialog_row_del(Thandle,Jhandle)
        rows = Jhandle.getComponent(0).getComponent(0).getSelectedRows+1;
        Ttmp=get(Thandle);
        if isempty(rows)
            rows = size(Ttmp.Data,1);
        end
        if rows > 1
            datatmp = Ttmp.Data(1:rows-1,:);
        end
        if rows < size(Ttmp.Data,1)
            datatmp(rows:size(Ttmp.Data,1)-1,:) = Ttmp.Data(rows+1:end,:);
        end
        if ~exist('datatmp','var')
            for nCol=1:length(Ttmp.ColumnEditable);
                if ~isempty(Ttmp.Data) & isnumeric(Ttmp.Data);
                    Ttmp.Data(1,nCol)=0;
                elseif strcmp(Ttmp.ColumnFormat{nCol},'numeric');
                    Ttmp.Data{1,nCol}=[];
                elseif strcmp(Ttmp.ColumnFormat{nCol},'logical');
                    Ttmp.Data{1,nCol}=false;
                else
                    Ttmp.Data{1,nCol}=' ';
                end;
            end;
        else
            Ttmp.Data = datatmp;
        end
        clearvars datatmp
        set(Thandle,'Data',Ttmp.Data);
        Ttmp = get(Thandle);
        Ttmp2 = get(gcf,'Position');
        Ttmp2(3) = Ttmp2(3)*Ttmp.Extent(3)/Ttmp.Position(3);
        Ttmp2(4) = Ttmp2(4)*Ttmp.Extent(4)/Ttmp.Position(4);
        if Ttmp2(3) > 1500
            Ttmp2(3) = 1500;
            Ttmp2(4) = Ttmp2(4)+22;
        end
        if Ttmp2(4) < 150
            Ttmp2(4) = 150;
        elseif Ttmp2(4) > 800
            Ttmp2(4) = 800;
            Ttmp2(2) = 70;
            Ttmp2(3) = Ttmp2(3)+22;
        end
        set(gcf,'Position',Ttmp2);
    end

    function generate_matrix(Thandle)
        matrix = get(Thandle,'Data');
        matrix = lab_generate_matrix(size(matrix,2));
        set(Thandle,'Data',matrix);
    end

    function read_montage(Thandle)
        [~,table] = lab_read_montage;
        if ~isempty(table)
            set(Thandle,'Data',table)
        end;
    end
% end support scripts table_dialog
% --------------------------------

end