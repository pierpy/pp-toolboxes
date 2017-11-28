% Helper script for TAPEEG to set 'FixedWidthFontName'
% (needed for linux systems were 'FixedWidthFontName' is sometimes wrong/empty)
%
% Written by F. Hatz 2014 University Hospital Basel

function lab_selectfont_fixedwidth(Force)

global Main_Path

if ~exist('Force','var')
    Force = false;
end

settings = [];
if ~exist(fullfile(Main_Path,'.ignore'),'dir')
    mkdir(fullfile(Main_Path,'.ignore'))
    if ispc
        fileattrib(fullfile(Main_Path,'.ignore'),'+h');
    end
end
File = fullfile(fullfile(Main_Path,'.ignore'),'FixedwidthFont');

if exist(File,'file')
    fid=fopen(File,'r');
    settings.font = fgetl(fid);
    fclose(fid);
    if ~isempty(find(strcmp(listfonts,settings.font),1)) & Force == false
        set(0,'FixedWidthFontName',settings.font);
        return
    else
        if isempty(find(strcmp(listfonts,settings.font),1))
            settings.font = get(0,'FixedWidthFontName');
        end
        if isempty(find(strcmp(listfonts,settings.font),1))
            settings.font = '';
        end
    end
else
    settings.font = get(0,'FixedWidthFontName');
    if isempty(find(strcmp(listfonts,settings.font),1))
        settings.font = '';
    end
end

% select font
Prompt = {'First run: Please select standard font for figures','';'Font','font'};
Formats.type = 'text';
Formats(2,1).type = 'none';
Formats(3,1).type = 'list';
Formats(3,1).style = 'popupmenu';
Formats(3,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(3,1).items = listfonts;
[settings,Cancelled] = inputsdlg(Prompt,'Select Font',Formats,settings);
while Cancelled == 1
    [settings,Cancelled] = inputsdlg(Prompt,'Select Font',Formats,settings);
end

% store font
fid = fopen(File,'w');
fprintf(fid,[settings.font '\n']);
fclose(fid);

% set font
set(0,'FixedWidthFontName',settings.font);