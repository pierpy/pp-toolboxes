function [cfg,skipprocessing] = lab_set_edit_markers(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header= [];
end
if ~exist('cfg','var') | ~isfield(cfg,'MARK')
    cfg.MARK = [];
end

% Convert old settings
if isfield(cfg.MARK,'markers') & ~isfield(cfg.MARK,'edit')
    cfg.MARK.edit = cfg.MARK.markers;
    cfg.MARK.markers = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Edit Markers','edit'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@edit_markers,'@ALL','@ALL'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Shuffle markers','shuffle'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'', 'markerinclude'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 220;
if exist('header','var') & isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
    tmp = unique(header.events.TYP);
    Formats(end,1).items = [{'all'} tmp(:)'];
else
    Formats(end,1).items = {'all'};
end

Prompt(end+1,:) = {'Rearrange','rearrange'};
Formats(end+1,1).type = 'check';

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Create markers','create'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'', 'markers'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 220;
if exist('header','var') & isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
    tmp = unique(header.events.TYP);
    Formats(end,1).items = [{'all'} tmp(:)'];
else
    Formats(end,1).items = {'all'};
end

Prompt(end+1,:) = {'Duration','duration'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999];
Formats(end,1).size = 80;

[cfg.MARK,Cancelled] = inputsdlg(Prompt,'Edit Markers',Formats,cfg.MARK);
if Cancelled == 1
    cfg.MARK = [];
    skipprocessing = 1;
end

  function settings = edit_markers(settings)
    if isempty(settings) | ~isfield(settings,'edit') | isempty(settings.edit)
        settings.edit = {'',0,'seconds',0,'add',false,'','all',false};
    elseif size(settings.edit,2) == 8
        settings.edit = [settings.edit(:,1:7) repmat({'all'},size(settings.edit,1),1) settings.edit(:,8)];
    end
    if exist('header','var') & isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
        markers = unique(header.events.TYP);
        markers = [markers(:)' cellstr('New Marker')];
    else
        markers = cellstr('New Marker');
    end
    Theader = {'Marker','Shift start','Unit','Change duration','Mode','New Marker','New Name','Selection','Delete'};
    ColumnFormat = {markers,'numeric',{'seconds','timeframes'},'numeric',{'add','fixed'},'logical','char',{'all','first','last'},'logical'};
    settings.edit = lab_table_dialog(settings.edit,Theader,'Edit Markers',2,ColumnFormat);
  end

end