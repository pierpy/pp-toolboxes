function [cfg,skipprocessing] = lab_set_plot_elec(cfg)

skipprocessing = 0;
global Main_Path

if ~exist('cfg','var')
    cfg = [];
end
if exist('cfg','var') & isstruct(cfg) & isfield(cfg,'x')
    tmp.LOCS = cfg;
    cfg = tmp;
    cfg.locs_file = '';
    clearvars tmp
end

% Read electrodes file
if ~isfield(cfg,'LOCS')
    if ~isfield(cfg,'loc_file') | ~exist(cfg.loc_file,'file')
        if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes.sfp'),'file')
            cfg.loc_file = fullfile(Main_Path,'electrodes.sfp');
        elseif ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes.els'),'file')
            cfg.loc_file = fullfile(Main_Path,'electrodes.els');
        elseif ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes.xyz'),'file')
            cfg.loc_file = fullfile(Main_Path,'electrodes.xyz');
        elseif exist(fullfile(pwd,'electrodes.sfp'),'file')
            cfg.loc_file = fullfile(pwd,'electrodes.sfp');
        elseif exist(fullfile(pwd,'electrodes.els'),'file')
            cfg.loc_file = fullfile(pwd,'electrodes.els');
        elseif exist(fullfile(pwd,'electrodes.xyz'),'file')
            cfg.loc_file = fullfile(pwd,'electrodes.xyz');
        else
            cfg.loc_file = '';
        end
    end
    if ~isempty(cfg.loc_file) & exist(cfg.loc_file,'file')
        cfg.LOCS = lab_read_locs(cfg.loc_file);
    else
        cfg.LOCS = [];
    end
end
if ~isfield(cfg,'loc_file')
    cfg.loc_file = '';
end

if ~isfield(cfg,'Mappings')
    MappingsList = lab_search(pwd,{'Mappings*.xls';'Mappings*.xlsx'},1);
    if ~isempty(MappingsList)
        cfg.Mappings = lab_read_mappings(MappingsList{1});
        if ~isempty(cfg.Mappings)
            cfg.Mappings.shortnames = true;
        end
    else
        cfg.Mappings = [];
    end
    clearvars MappingsList
end

if ~isempty(cfg.LOCS) & ~isempty(cfg.Mappings) & isfield(cfg.Mappings,'mappingsChannelsFile') & ...
        isfield(cfg.LOCS,'x') & cfg.Mappings.mappingsChannelsFile ~= length(cfg.LOCS.x)
    excludeM = lab_get_exclude(cfg.Mappings.mappingsChannelsFile);
    excludeL = lab_get_exclude(length(cfg.LOCS.x));
    if ~isempty(excludeM) & length(setdiff(1:cfg.Mappings.mappingsChannelsFile,excludeM)) == length(cfg.LOCS.x)
        cfg.Mappings = lab_reduce_mappings(cfg.Mappings,exclude);
    elseif ~isempty(excludeL) & length(setdiff(1:length(cfg.LOCS.x),excludeM)) == cfg.Mappings.mappingsChannelsFile
        cfg.LOCS = lab_reduce_locs(cfg.LOCS,excludeL);
    else
        cfg.Mappings = [];
    end
end

if ~isfield(cfg,'plothead')
    cfg.plothead = true;
end
if ~isfield(cfg,'onlyhead')
    cfg.plothead = false;
end

List_locs = {};
if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes'),'dir')
    List_locs = lab_search(fullfile(Main_Path,'electrodes'),{'*.els','*.xyz','*.sfp','*.elc','*.spi'},true,true,1);
end
List_locs = cat(1,{'Select File'},List_locs(:));
List_locs = cat(1,{cfg.loc_file},List_locs(:));

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'LOC-file','loc_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = List_locs;
Formats(end,1).size = 250;
Formats(end,1).callback = {@load_locs,'@ALL','@ALL'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'LOCS','LOCS'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [160 160];
Formats(end,1).enable = 'inactive';
Formats(end,1).span = [8 2];

Prompt(end+1,:) = {'Mappings','Mappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_mappings,'@ALL','@ALL'};

Prompt(end+1,:) = {'Convert',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@convert_mappings,'@ALL','@ALL'};

Prompt(end+1,:) = {'Match electrodes & input-data','Match'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@lab_plot_match_numchans,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Plot head','plothead'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Restrict to head','onlyhead'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'  ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [5 1];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Data','DATA'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@read_data,'@ALL','@ALL'};

Prompt(end+1,:) = {'settings','PLOT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_settings,'@ALL','@ALL'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Save Picture','Store'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'','DATA_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).limits = [1 0];
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Plot',''};
Formats(end+1,1).type = 'button';
Formats(end,1).size = 90;
Formats(end,1).callback = {@lab_plot_elec_plotall,'','@ALL'};

Options.CancelButton = 'off';
Options.ButtonNames = {'Close'};
[cfg,Cancelled] = inputsdlg(Prompt,'PLOT Locs',Formats,cfg,Options);
if Cancelled == 1
    cfg = [];
    skipprocessing = 1;
    return
end

    function settings = load_locs(settings)
        if ~isempty(settings.loc_file) & strcmp(settings.loc_file,'Select File')
            [LOC_file,LOC_filepath] = uigetfile('*.els;*.sfp;*.xyz','Select LOC-file');
            if LOC_file == 0
                return
            end
            settings.loc_file = fullfile(LOC_filepath,LOC_file);
        end
        if exist(settings.loc_file,'file')
            settings.LOCS = lab_read_locs(settings.loc_file);
            if ~isempty(settings.LOCS)
                if min(settings.LOCS.z) == 0 & max(settings.LOCS.z) == 0
                    settings.plothead = false;
                else
                    settings.plothead = true;
                end
                settings = lab_plot_check_numchans(settings);
                if settings.Match == false
                    settings = lab_plot_match_numchans(settings);
                end
            else
                settings.plothead = false;
                settings.Match = false;
            end
        else
            settings.LOCS = [];
            settings.plothead = false;
            settings.Match = false;
        end
    end
    
    function settings = convert_mappings(settings)
        if ~isempty(settings.LOCS) & ~isempty(settings.Mappings)
            settings.LOCS = lab_Mappings2Locs(settings.Mappings,settings.LOCS);
            settings.Mappings = [];
        end
    end
    
    function settings = read_data(settings)
        if isfield(settings,'Mappings') & isfield(settings.Mappings,'mappingstitle')
            labels = settings.Mappings.mappingstitle;
        elseif isfield(settings,'LOCS') & isfield(settings.LOCS,'labels')
            labels = settings.LOCS.labels';
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
        [DATA,settings] = lab_plot_read_data(settings.DATA,settings,labels,Mode);
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
        settings.PLOT = lab_plot_settings(settings.PLOT,settings.DATA,0);
    end
    
    function settings = load_mappings(settings)
        settings.Mappings = lab_load_mappings(settings.Mappings,settings,'Mappings.xls',settings.LOCS,true);
        if ~isempty(settings.Mappings) & ~isempty(settings.LOCS) & settings.Mappings.mappingsChannels ~= size(settings.LOCS.x,2)
            settings.Mappings = [];
        end
        settings = lab_plot_check_numchans(settings);
    end
end