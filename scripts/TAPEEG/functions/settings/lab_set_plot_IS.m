function [settings,skipprocessing] = lab_set_plot_IS(settings)

skipprocessing = 0;
global Main_Path

% Preload brain template
if exist(fullfile(pwd,'MNIbrain.mat'),'file')
    disp(['    load ' fullfile(pwd,'MNIbrain.mat')])
    settings.Brain = load(fullfile(pwd,'MNIbrain.mat'));
    settings.Brain.Brain_file = fullfile(pwd,'MNIbrain.mat');
    settings.Labels = char(settings.Brain.labels);
elseif ~isempty(Main_Path) & exist(fullfile(Main_Path,'MNIbrain.mat'),'file')
    disp(['    load ' fullfile(Main_Path,'MNIbrain.mat')])
    settings.Brain = load(fullfile(Main_Path,'MNIbrain.mat'));
    settings.Brain.Brain_file = fullfile(Main_Path,'MNIbrain.mat');
    settings.Labels = char(settings.Brain.labels);
end

% Read electrodes file

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Template-Regions','Labels'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).limits = [0 8];
Formats(end,1).size = [170 180];
Formats(end,1).enable = 'inactive';
Formats(end,1).span = [8 2];

Prompt(end+1,:) = {'Brain-Template','Brain'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@show_brain,'@ALL','@ALL'};

Prompt(end+1,:) = {'Load Brain',''};
Formats(end+1,1).type = 'button';
Formats(end,1).size = 90;
Formats(end,1).callback = {@load_brain,'@ALL','@ALL'};

Prompt(end+1,:) = {'Generate Brain',''};
Formats(end+1,1).type = 'button';
Formats(end,1).size = 90;
Formats(end,1).callback = {@create_brain,'@ALL','@ALL'};

Prompt(end+1,:) = {'Select regions',''};
Formats(end+1,1).type = 'button';
Formats(end,1).size = 90;
Formats(end,1).callback = {@select_regions,'@ALL','@ALL'};

Prompt(end+1,:) = {'Mappings',''};
Formats(end+1,1).type = 'button';
Formats(end,1).size = 90;
Formats(end,1).callback = {@load_mappings,'@ALL','@ALL'};

Prompt(end+1,:) = {'Match regions & input-data','Match'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@lab_plot_check_numchans,'@ALL','@ALL'};

Prompt(end+1,:) = {'  ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [3 1];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Data','DATA'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@read_data,'@ALL','@ALL'};

Prompt(end+1,:) = {'settings','PLOT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_settings,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Save Picture','Store'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'','DATA_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).limits = [1 0];
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Plot',''};
Formats(end+1,1).type = 'button';
Formats(end,1).size = 90;
Formats(end,1).callback = {@lab_plot_IS_plotall,'','@ALL'};

Options.CancelButton = 'off';
Options.ButtonNames = {'Close'};
[settings,Cancelled] = inputsdlg(Prompt,'PLOT IS',Formats,settings,Options);
if Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
end
    
end

function show_brain(settings)
    Brain = settings.Brain;
    Labels = Brain.labels;
    Color = rand(length(Labels),3);
    Color(end+1,:) = 1;
    Facecolor = Brain.mapsall * Color;
    F = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none','Name','Bain-Template');
    m1 = uimenu(F,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(F,'Label','Edit');
    uimenu(m2,'Label','Rotate on/off','Callback',@(~,~)rotate3d);
    lab_show_brain(F,{Brain.faces},{Brain.vertices},{Facecolor},lines,'T',0,1);
    light('Position',[0 0 -1]);
    for i = 1:length(Labels)
        Legend(i).Color = Color(i,:); %#ok<AGROW>
        Legend(i).Text = regexprep(Labels{i},'_',' '); %#ok<AGROW>
        Legend(i).Mode = 'Legend'; %#ok<AGROW>
    end
    lab_plot_legend(Legend);
    rotate3d on
end

function settings = load_brain(settings)
    [Brain_file,Brain_filepath] = uigetfile('*.mat','Select Brain-Template (MNIbrain.mat)');
    if Brain_file == 0
        settings.Brain = [];
        settings.Match = false;
        settings.Labels = '';
        return
    end
    Brain_file = fullfile(Brain_filepath,Brain_file);
    if exist(Brain_file,'file') & strcmp(Brain_file(end-2:end),'mat')
        settings.Brain = load(Brain_file);
        settings.Brain.Brain_file = Brain_file;
        settings.Match = false;
        settings.Labels = char(settings.Brain.labels);
    else
        settings.Brain = [];
        settings.Match = false;
        settings.Labels = '';
    end
end

function settings = create_brain(settings)
    settings.Brain = lab_calculate_showbrain('MRI-file');
    if isempty(settings.Brain)
        settings.Labels = '';
        settings.Match = false;
    else
        settings.Labels = char(settings.Brain.labels);
        settings.Match = false;
    end
end

function settings = select_regions(settings)
    Brain = lab_calculate_showbrain(settings.Brain,1);
    if ~isempty(Brain)
        settings.Brain = Brain;
        settings.Labels = char(settings.Brain.labels);
        settings.Match = false;
    end
end

function settings = load_mappings(settings)
    Brain = lab_brain_template2regions(settings.Brain);
    if ~isempty(Brain)
        settings.Brain = Brain;
        settings.Labels = char(settings.Brain.labels);
        settings.Match = false;
    end
end

function settings = read_data(settings)
    if isfield(settings,'Brain') & isfield(settings.Brain,'labels')
        labels = settings.Brain.labels(:);
    else
        labels = [];
    end
    if isempty(settings.DATA)
        Mode = 'Select';
    elseif isfield(settings.DATA,'labels') & settings.DATA(1).labels == true
        Mode = 'Labels';
    else
        Mode = 'Add';
    end
    if isfield(settings.Brain,'regions') & ~isempty(settings.Brain.regions)
        doIS = 2;
    else
        doIS = 1;
    end
    [DATA,settings] = lab_plot_read_data(settings.DATA,settings,labels,Mode,doIS);
    settings.DATA = DATA;
    if isempty(DATA)
        settings.PLOT = [];
        settings.Match = false;
        return
    else
        settings = lab_plot_check_numchans(settings);
    end
end

function settings = set_settings(settings)
    settings = lab_plot_check_numchans(settings);
    if isfield(settings.Brain,'regions') & ~isempty(settings.Brain.regions)
        doIS = 2;
    else
        doIS = 1;
    end
    settings.PLOT = lab_plot_settings(settings.PLOT,settings.DATA,doIS);
end