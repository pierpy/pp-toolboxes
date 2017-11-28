function mappings = lab_reduce_mappings(mappings,exclude,cfg)

if ~exist('cfg','var')
    cfg = [];
end

selection = 1:mappings.mappingsChannelsFile;
if isfield(mappings,'mappingsexclude')
    selection = setdiff(selection,mappings.mappingsexclude);
    excludeB = mappings.mappingsexclude;
else
    excludeB = [];
end

mappingstmp = zeros(size(mappings.mappings,2),mappings.mappingsChannels);
for i = 1:size(mappings.mappings,2)
    mappingstmp(i,mappings.mappings{1,i}) = 1;
end

if ~exist('exclude','var') | isempty(exclude)
    if isfield(cfg,'exclude')
        exclude = cfg.exclude;
        [~,exclude] = intersect(selection,exclude);
    else
        if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
            if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
                exclude = cfg.EXTRA.exclude;
                [~,exclude] = intersect(selection,exclude);
            else
                exclude = [];
            end
        else
            strlist = cellstr(num2str(selection'));
            if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
                strdefault = cfg.EXTRA.exclude;
                [~,strdefault] = intersect(selection,strdefault);
            elseif ~isempty(lab_get_exclude(mappings.mappingsChannelsFile))
                strdefault = lab_get_exclude(mappings.mappingsChannelsFile);
                [~,strdefault] = intersect(selection,strdefault);
            else
                strdefault = [];
            end
            exclude = listdlg('PromptString','Excluded channels (mappings-file):','SelectionMode','multiple', ...
                'ListString',strlist,'InitialValue',strdefault,'CancelString','None','ListSize',[220 350]);
            clearvars strlist strdefault
        end
    end
end

if ~isempty(exclude) & exclude ~= 0
    includechans = setdiff((1:size(mappingstmp,2)),exclude);
    mappingstmp = mappingstmp(:,includechans);
    clearvars includechans
    for i = 1:size(mappings.mappings,2)
        mappings.mappings{i} = find(mappingstmp(i,:) == 1);
    end
    mappings.mappingsChannels = size(mappingstmp,2);
    mappings.mappingsexclude = union(selection(exclude),excludeB);
    mappings.mappingsexclude = mappings.mappingsexclude(:)';
end