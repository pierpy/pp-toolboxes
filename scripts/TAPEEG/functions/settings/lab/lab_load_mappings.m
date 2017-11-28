% Read xls files with channel mappings info
%
% Mappings = lab_read_mappings(cfg,ISflag)
%
% cfg    = structure with config (optional)
% ISflag = as default the file 'MappingsIS.xls' is loaded
%          (in place of 'Mappings.xls')
%
% written by F. Hatz 2012

function [Mappings,Mappings_file] = lab_load_mappings(Mappings,cfg,filename,Locs,flagnames)

global Main_Path
    
if ~exist('flagnames','var')
    flagnames = false;
end
if ~exist('Locs','var')
    Locs = [];
end
if ~exist('filename','var') | isempty(filename)
    filename = 'Mappings.xls';
end

if ~exist('cfg','var') | ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end
if ~exist('Mappings','var')
    Mappings = [];
end

if isempty(Mappings)
    if exist(fullfile(cfg.settings_path,filename),'file')
        mappings_file = fullfile(cfg.settings_path,filename);
    elseif exist(fullfile(cfg.settings_path,[filename 'x']),'file')
        mappings_file = fullfile(cfg.settings_path,[filename 'x']);
    elseif exist(fullfile(pwd,filename),'file')
        mappings_file = fullfile(pwd,filename);
    elseif exist(fullfile(pwd,[filename 'x']),'file')
        mappings_file = fullfile(pwd,[filename 'x']);
    end
    if exist('mappings_file','var')
        Mappings = read_mappings(mappings_file,cfg);
    end
end

if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
    settings.mappings = Mappings;
    if flagnames == true & isfield(Mappings,'shortnames')
        settings.shortnames = Mappings.shortnames;        
    end
    if flagnames == true & isfield(Mappings,'plotdots')
        settings.plotdots = Mappings.plotdots;        
    end
    
    if exist(fullfile(Main_Path,'mapping'),'dir')
        List_mappings = lab_search(fullfile(Main_Path,'mapping'),{'*.xls','*.xlsx'},true,true,1);
        List_mappings = cat(1,{'Select File'},List_mappings(:));
    else
        List_mappings = {'Select File'};
    end
    
    if exist('mappings_file','var')
        settings.mappings_file = mappings_file;
        List_mappings = cat(1,{mappings_file},List_mappings(:));
    else
        settings.mappings_file = '';
        List_mappings = cat(1,{''},List_mappings(:));
    end
    
    Prompt = cell(0,2);
    Formats = [];
    
    if ~exist('ISflag','var')
        Prompt(end+1,:) = {'Mappings-File','mappings_file'};
    else
        Prompt(end+1,:) = {'Mappings-IS-File','mappings_file'};
    end
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = List_mappings;
    Formats(end,1).callback =  {@read_mappings,{'mappings','mappings_file'},'mappings_file',cfg};
    Formats(end,1).size = 270;
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Mappings','mappings'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'result';
    Formats(end,1).size = [190 90];
    Formats(end,1).enable = 'inactive';
    Formats(end,1).callback = {@edit_mappings,'mappings','mappings',Locs};
    Formats(end,1).span = [2 1];
    
    Prompt(end+1,:) = {'Reduce',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@lab_reduce_mappings,'mappings','mappings',[],cfg};
    
    Prompt(end+1,:) = {'Edit',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@edit_mappings,'mappings','mappings',Locs};
    
    if flagnames == true
        Prompt(end+1,:) = {'Use short names','shortnames'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).span = [1 2];
        
        Prompt(end+1,:) = {'Plot as dots','plotdots'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).span = [1 2];
    end
    
    [settings,Cancelled] = inputsdlg(Prompt,'Load Mappings',Formats,settings);
    if Cancelled == 1
        Mappings = [];
        Mappings_file = '';
        return
    end
    Mappings = settings.mappings;
    if flagnames == true
        Mappings.shortnames = settings.shortnames;
        Mappings.plotdots = settings.plotdots;
    end
    Mappings_file = settings.mappings_file;
end

end

function [Mappings,Mappings_file] = read_mappings(Mappings_file,cfg)
    if ~isempty(Mappings_file) & strcmp(Mappings_file,'Select File')
        [Mappings_file,Mappings_path] = uigetfile({'*.xls;*.xlsx','Excel-File'},'Select Mappings-File');
        Mappings_file = fullfile(Mappings_path,Mappings_file);
    end
    if isempty(Mappings_file) | ~exist(Mappings_file,'file')
        Mappings_file = '';
        Mappings = [];
        return
    end
    Mappings = lab_read_mappings(Mappings_file);
    if isempty(Mappings)
        return
    end
    if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans') & ...
            ~isempty(cfg.EXTRA.numdatachans) & cfg.EXTRA.numdatachans ~= Mappings.mappingsChannels
        if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
            exclude = cfg.exclude;
        else
            exclude = lab_get_exclude(Mappings.mappingsChannelsFile);
        end
        if length(setdiff(1:Mappings.mappingsChannelsFile,exclude)) == cfg.EXTRA.numdatachans
            Mappings = lab_reduce_mappings(Mappings,exclude);
        end
    end
end

function MappingsOut = edit_mappings(Mappings,Locs)
  if isempty(Mappings) | ~isfield(Mappings,'mappings')
      MappingsOut = [];
      return
  end
  if ~exist('Locs','var') | ~isfield(Locs,'x')
      MapEdit = [];
      for i = 1:size(Mappings.mappings,2)
          MapEdit(i).Name = Mappings.mappingstitle{i,1};
          MapEdit(i).Channels = Mappings.mappings{1,i};
      end
      FORMATS{strcmp(fieldnames(MapEdit),'Name')} = 'string';
      FORMATS{strcmp(fieldnames(MapEdit),'Channels')} = cellstr(num2str((1:Mappings.mappingsChannels)'));
      
      MapEdit = inputsdlg(MapEdit,'Edit mappings',FORMATS);
      
      MappingsOut.mappingsChannels = Mappings.mappingsChannels;
      MappingsOut.mappingsChannelsFile = Mappings.mappingsChannelsFile;
      for i = 1:length(MapEdit)
          MappingsOut.mappings{1,i} = MapEdit(i).Channels;
          MappingsOut.mappingstitle{i,1} = MapEdit(i).Name;
          tmp = MappingsOut.mappingstitle{i,1};
          tmp=textscan(tmp,'%s');
          tmp = tmp{1,1};
          title = '';
          for j = 1:size(tmp,1);
              if size(tmp,1) == 1 & length(tmp{1,1}) > 1
                  title = [upper(tmp{1,1}(1)) lower(tmp{1,1}(2))];
              elseif j < 4
                  title = [title upper(tmp{j,1}(1))];
              end
          end
          MappingsOut.mappingstitleS{i,1} = title;
      end
  else
      if Mappings.mappingsChannels > size(Locs.x,2)
          MappingsOut = Mappings;
          disp('   plotting of mappings not possible, mismatch mappings & locs')
          return
      end
      plot.indexed = Mappings;
      plot.Color = [1 1 1];
      plot.ColorIdx = [1 0 0];
      plot.LOCS = Locs;
      plot.Title = 'Mappings';
      [indexed,mappingstmp] = lab_plot_locs(plot,1,0,1);
      if ~isempty(indexed)
          MappingsOut.mappingsChannelsFile = size(Locs.x,2);
          MappingsOut.mappingsChannels = size(Locs.x,2);
          MappingsOut.mappings = indexed;
          for i = 1:size(indexed,2)
              if ~isempty(mappingstmp)
                  MappingsOut.mappingstitle{i,1} = mappingstmp(1,i).Name;
              else
                  MappingsOut.mappingstitle{i,1} = ['Mapping' num2str(i)];
              end
          end
          for i = 1:size(MappingsOut.mappingstitle,1)
              tmp = MappingsOut.mappingstitle{i,1};
              tmp=textscan(tmp,'%s');
              tmp = tmp{1,1};
              title = '';
              for j = 1:size(tmp,1);
                  if size(tmp,1) == 1 & length(tmp{1,1}) > 1
                      title = [upper(tmp{1,1}(1)) lower(tmp{1,1}(2))];
                  elseif j < 4
                      title = [title upper(tmp{j,1}(1))];
                  end
              end
              MappingsOut.mappingstitleS{i,1} = title;
          end
      else
          MappingsOut = [];
      end
  end
end