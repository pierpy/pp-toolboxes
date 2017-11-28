function settings = lab_set_plot_locs(settings)

global Main_Path

if ~exist('settings','var') | ~isfield(settings,'color')
    if exist(fullfile(Main_Path,'electrodes.els'),'file')
        [settings.LOCS,settings.LOC_file] = load_locs(fullfile(Main_Path,'electrodes.els'));
    elseif exist(fullfile(Main_Path,'electrodes.sfp'),'file')
        [settings.LOCS,settings.LOC_file] = load_locs(fullfile(Main_Path,'electrodes.sfp'));
    elseif exist(fullfile(Main_Path,'electrodes.xyz'),'file')
        [settings.LOCS,settings.LOC_file] = load_locs(fullfile(Main_Path,'electrodes.xyz'));
    else
        settings.LOCS = [];
        settings.LOC_file = '';
    end
    settings.color = [1 1 1];
    settings.indexed = [];
    settings.colori = [1 0 0];
    settings.colormode = 'color';
end

List_locs = {};
if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes'),'dir')
    List_locs = lab_search(fullfile(Main_Path,'electrodes'),{'*.els','*.xyz','*.sfp','*.elc','*.spi'},true,true,1);
end
List_locs = cat(1,{'Select File'},List_locs(:));
List_locs = cat(1,{settings.LOC_file},List_locs(:));

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'LOC-file','LOC_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = List_locs;
Formats(end,1).size = 250;
Formats(end,1).callback = {@load_locs,{'LOCS','LOC_file'},'LOC_file'};

Prompt(end+1,:) = {'LOCS','LOCS'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [250 160];
Formats(end,1).enable = 'inactive';

Prompt(end+1,:) = {'Color','color'};
Formats(end+1,1).type = 'color';

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Indexed','indexed'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 250;
Formats(end,1).limits = [-inf inf]; % if not [0 1] data = numeric

Prompt(end+1,:) = {'Color (indexed)','colori'};
Formats(end+1,1).type = 'color';

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Colormode','colormode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'color','gray'};
Formats(end,1).callback = {@set_colormode,'@ALL','@ALL'};

[settings,Cancelled] = inputsdlg(Prompt,'PLOT Locs',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    return
end

    function [LOCS,LOC_file] = load_locs(LOC_file)
        if ~exist(LOC_file,'file') | strcmp(LOC_file,'Select File')
            [LOC_file,LOC_filepath] = uigetfile('*.els;*.sfp;*.xyz','Select LOC-file');
            if LOC_file == 0
                return
            end
            LOC_file = fullfile(LOC_filepath,LOC_file);
        end
        if ~isempty(LOC_file) & exist(LOC_file,'file')
            LOCS = lab_read_locs(LOC_file);
        else
            LOCS = [];
        end
    end
    function settings = set_colormode(settings)
        if strcmp(settings.colormode,'color')
            settings.color = [1 1 1];
            settings.colori = [1 0 0];
        else
            settings.color = [1 1 1];
            settings.colori = [0.5 0.5 0.5];
        end
    end
end
