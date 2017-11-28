% Load LOCS-file for processing
%
% LOCS = lab_load_locs(LOCS,cfg.numchans)
%
% cfg:  structure with config (optional)
%
% written by F. Hatz 2012

function LOCS = lab_load_locs(LOCS,cfg,numchans)

global Main_Path

filename = 'electrodes.els';

if ~exist('cfg','var') | ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end
if ~exist('LOCS','var')
    LOCS = [];
end
if ~exist('numchans','var') | isempty(numchans) | numchans == 0
    if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & ...
            ~isempty(cfg.EXTRA.numdatachans)
        numchans = cfg.EXTRA.numdatachans;
    else
        numchans = [];
    end
end
      
if isempty(LOCS)
    if exist(fullfile(cfg.settings_path,filename),'file')
        LOCS_file = fullfile(cfg.settings_path,filename);
    elseif exist(fullfile(pwd,filename),'file')
        LOCS_file = fullfile(pwd,filename);
    end
    if exist('LOCS_file','var')
        LOCS = lab_read_locs(LOCS_file);
    end
end

if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
    List_locs = {};
    if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes'),'dir')
        List_locs = lab_search(fullfile(Main_Path,'electrodes'),{'*.els','*.xyz','*.sfp','*.elc','*.spi'},true,true,1);
    end
    List_locs = cat(1,{'Select File'},List_locs(:));

    settings.LOCS = LOCS;
    if exist('LOCS_file','var')
        settings.LOCS_file = LOCS_file;
        List_locs = cat(1,{LOCS_file},List_locs(:));
    else
        settings.LOCS_file = [];
        List_locs = cat(1,{''},List_locs(:));
    end
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'LOCS-File','LOCS_file'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = List_locs;
    Formats(end,1).callback =  {@read_locs,{'LOCS','LOCS_file'},'LOCS_file',cfg,numchans};
    Formats(end,1).size = 300;
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'LOCS','LOCS'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'result';
    Formats(end,1).size = [240 110];
    Formats(end,1).enable = 'inactive';
    Formats(end,1).span = [2 1];
    Formats(end,1).callback = {@edit_locs,'LOCS','LOCS'};
    
    Prompt(end+1,:) = {'Reduce',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@lab_reduce_locs,'LOCS','LOCS',[],cfg};

    Prompt(end+1,:) = {'Edit',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@edit_locs,'LOCS','LOCS'};
    
    [settings,Cancelled] = inputsdlg(Prompt,'Load LOCS',Formats,settings);
    if Cancelled == 1
        LOCS = [];
        return
    end
    LOCS = settings.LOCS;
end

end

function [LOCS,LOCS_file] = read_locs(LOCS_file,cfg,numchans)
    if ~isempty(LOCS_file) & strcmp(LOCS_file,'Select File')
        [LOCS_file,LOCS_path] = uigetfile({'*.els;*.xyz;*.sfp;*.elc;*.spi','LOCS-File'},'Select LOCS-File');
        LOCS_file = fullfile(LOCS_path,LOCS_file);
    end
    if isempty(LOCS_file) | ~exist(LOCS_file,'file')
        LOCS_file = '';
        LOCS = [];
        return
    end
    LOCS = lab_read_locs(LOCS_file);
    if isempty(LOCS)
        return
    end
    if ~isempty(numchans) & numchans ~= length(LOCS.x)
        if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
            exclude = cfg.exclude;
        else
            exclude = lab_get_exclude(length(LOCS.x));
        end
        if length(setdiff(1:length(LOCS.x),exclude)) == numchans
            LOCS = lab_reduce_locs(LOCS,exclude);
        end
    end
end

function Locs = edit_locs(Locs)
  if ~exist('Locs','var') | isempty(Locs)
      return
  end
  plot.indexed = [];
  plot.Color = [1 1 1];
  plot.ColorIdx = [1 0 0];
  plot.LOCS = Locs;
  plot.Title = 'Locs';
  [~,~,Locs] = lab_plot_locs(plot,1,1,0);
end
