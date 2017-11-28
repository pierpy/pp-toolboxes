function settings = lab_get_eegsource(settings,cfg,header,Filename)

if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('settings','var')
    settings = [];
end

if any(strcmp(settings.eegsource,'laplacian'))
    if (~isfield(settings,'LAPL') | isempty(settings.LAPL))
        [settings.LAPL,skipprocessing] = lab_get_laplacian(settings.LAPL);
        if skipprocessing == 1
            settings.LAPL = [];
            if ischar(settings.eegsource)
                settings.eegsource = 'mean';
            else
                tmp = ~strcmp(settings.eegsource,'laplacian');
                if max(tmp) == 1
                    settings.eegsource = settings.eegsource(tmp);
                else
                    settings.eegsource = {'mean'};
                end
            end
        end
    end
else
    settings.LAPL = [];
end
if any(strcmp(settings.eegsource,'montage'))
    if (~isfield(settings,'montage') | isempty(settings.montage))
        if ~isfield(settings,'montage')
            settings.montage = [];
        end
        if exist('Filename','var')
            settings.montage = lab_load_montage(settings.montage,cfg,header,Filename);
        else
            settings.montage = lab_load_montage(settings.montage,cfg,header);
        end
        if isempty(settings.montage)
            if ischar(settings.eegsource)
                settings.eegsource = 'mean';
            else
                tmp = ~strcmp(settings.eegsource,'montage');
                settings.eegsource = settings.eegsource(tmp);
            end
        end
    end
else
    settings.montage = [];
end
if any(strcmp(settings.eegsource,'channels'))
    if exist('header','var') & isfield(header,'channels')
        strlist = cellstr(header.channels);
        [selection] = listdlg('PromptString','Reference channels','SelectionMode','multiple', ...
            'ListString',strlist,'ListSize',[230 200],'CancelString','None');
        if ~isempty(selection)
            if ischar(settings.eegsource)
                settings.eegsource = num2str(selection);
            else
                tmp = strcmp(settings.eegsource,'channels');
                settings.eegsource{tmp} = num2str(selection);
                clearvars tmp
            end
        else
            if ischar(settings.eegsource)
                settings.eegsource = 'none';
            else
                tmp = ~strcmp(settings.eegsource,'channels');
                settings.eegsource = settings.eegsource(tmp);
                clearvars tmp
            end
        end
    else
        tmp.eegsource = '';
        Prompt = cell(0,2);
        Formats = [];
        Prompt(end+1,:) = {'Channels for reference','eegsource'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'text';
        Formats(end,1).size = 90; % automatically assign the height
        [tmp,Cancelled] = inputsdlg(Prompt,'Reference channels',Formats,tmp);
        if Cancelled == 1
            if ischar(settings.eegsource)
                settings.eegsource = 'mean';
            else
                tmp2 = ~strcmp(settings.eegsource,'channels');
                settings.eegsource = settings.eegsource(tmp2);
                clearvars tmp2
            end
        else
            if ischar(settings.eegsource)
                settings.eegsource = tmp.eegsource;
            else
                tmp2 = strcmp(settings.eegsource,'channels');
                settings.eegsource{tmp2} = tmp.eegsource;
                clearvars tmp2
            end
        end
        clearvars tmp
    end
end
if isfield(settings,'excludebad') & ~any(strcmp(settings.eegsource,'input'))
    settings.excludebad = false;
end
