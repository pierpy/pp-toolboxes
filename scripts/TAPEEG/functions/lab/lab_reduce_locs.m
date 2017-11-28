function LOCS = lab_reduce_locs(LOCS,exclude,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('LOCS','var') | isempty(LOCS) | ~isfield(LOCS,'x')
    return
end
if ~exist('exclude','var') | isempty(exclude)
    if isfield(cfg,'exclude')
        exclude = cfg.exclude;
    else
        if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
            if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
                exclude = cfg.EXTRA.exclude;
            else
                exclude = [];
            end
        else
            strlist = cellstr(num2str((1:length(LOCS.x))'));
            if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
                strdefault = cfg.EXTRA.exclude;
            elseif ~isempty(lab_get_exclude(length(LOCS.x)))
                strdefault = lab_get_exclude(length(LOCS.x));
            else
                strdefault = [];
            end
            exclude = listdlg('PromptString','Excluded channels (LOCS-file):','SelectionMode','multiple', ...
                'ListString',strlist,'InitialValue',strdefault,'CancelString','None','ListSize',[220 350]);
            clearvars strlist strdefault
        end
    end
end

if ~isempty(exclude)
    include = 1:length(LOCS.x);
    include = setdiff(include,exclude);
    LOCS.x = LOCS.x(1,include);
    LOCS.y = LOCS.y(1,include);
    LOCS.z = LOCS.z(1,include);
    LOCS.labels = LOCS.labels(1,include);
    LOCS = lab_locs2sph(LOCS);
end