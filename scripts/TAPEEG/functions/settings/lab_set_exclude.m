function cfg = lab_set_exclude(cfg,header)

global THeader

if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'exclude')
    cfg.exclude = [];
end

disp ('Ask for channels to exclude')
if isfield(header,'locs') & ~isempty(header.locs)
    exclude = lab_load_exclude(cfg.exclude,cfg,header,1);
    settings.indexed = exclude;
    settings.Color = [1 1 1];
    settings.ColorIdx = [1 0 0];
    settings.LOCS = header.locs;
    settings.Title = 'Excluded channels';
    cfg.exclude = lab_plot_locs(settings,1,0,0,0);
else
    cfg.exclude = lab_load_exclude(cfg.exclude,cfg,header);
end

if ~isempty(cfg.exclude) & isfield(header,'numdatachannels')
    cfg.exclude = cfg.exclude(cfg.exclude<=header.numdatachannels);
end

if ~isempty(cfg.exclude)
    cfg.EXTRA.exclude = cfg.exclude;
elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
    cfg.EXTRA = rmfield(cfg.EXTRA,'exclude');
end