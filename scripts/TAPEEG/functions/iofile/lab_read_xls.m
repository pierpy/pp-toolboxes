% result = lab_read_xls(filename)
% 
% by F. Hatz, Neurology Basel

function result = lab_read_xls(Filename,skipselect)

if ~exist('skipselect','var')
    skipselect = false;
end

if ~exist('Filename','var')
    [Filename,Filepath] = uigetfile('*.xls;*.xlsx','Select xls-file');
    Filename = fullfile(Filepath,Filename);
end

result = [];
if strfind(Filename,'Mappings') | strfind(Filename,'mappings')
    try
        result = lab_read_mappings(Filename);
    catch %#ok<CTCH>
        disp('   Invalid format for mappings-file')
    end
elseif strfind(Filename,'Montage') | strfind(Filename,'montage')
    try
        result = lab_read_montage(Filename);
    catch %#ok<CTCH>
        disp('   Invalid format for montage-file')
    end
end
if ~isempty(result)
    return
end

[~,sheets] = xlsfinfo(Filename);
if skipselect == true
    sheet = sheets{1};
    selection = 1;
elseif length(sheets) == 1
    sheet = sheets{1};
    selection = 1;
elseif length(sheets) == 3 & strcmp(sheets{1},'Tabelle1') & strcmp(sheets{2},'Tabelle2') & strcmp(sheets{3},'Tabelle3')
    sheet = sheets{1};
    selection = 1;
else
    selection = listdlg('PromptString','Select table','SelectionMode', ...
        'single','ListString',sheets,'ListSize',[150 170]);
    if isempty(selection)
        return
    else
        pause(0.2);
    end
    sheet = sheets{selection};
end
if ispc
    [~,~,result] = xlsread(Filename,sheet);
else
    [~,~,result] = xlsread(Filename,selection,'','basic');
end

if ~isempty(result) & iscell(result)
    tmp = cellfun(@testnan,result);
    tmp1 = find(sum(tmp,2) < size(tmp,2));
    tmp2 = find(sum(tmp,1) < size(tmp,1));
    if isempty(tmp1) | isempty(tmp2)
        result = [];
    else
        result = result(tmp1,tmp2);
    end
    result(cellfun(@isempty,result)) = {NaN};
end

if ~isempty(result) & iscell(result) & isnumeric(result{1,1})
    if size(result,2) > 1 & ~isnumeric(result{1,2})
        return
    elseif size(result,1) > 1 & ~isnumeric(result{2,1})
        return
    end
    try %#ok<TRYNC>
        result = cell2mat(result);
    end
end

end

function output = testnan(input)
    if isnumeric(input) & isnan(input)
        output = true;
    else
        output = false;
    end
end

