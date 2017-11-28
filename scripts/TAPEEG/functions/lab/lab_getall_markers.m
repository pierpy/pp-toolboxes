function Marker = lab_getall_markers(header,cfg)

if isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
    Marker = unique(header.events.TYP);
elseif isfield(cfg,'TRIGGER') & ~isempty(cfg.TRIGGER)
    Marker = cfg.TRIGGER;
else
    Marker = [];
end

if isfield(cfg,'MARK') & isfield(cfg.MARK,'edit') & ~isempty(cfg.MARK.edit) & ...
        size(cfg.MARK.edit,2) >= 7
    Marker = union(Marker,cfg.MARK.edit(:,7));
end

if isfield(cfg,'BADELEC') & isfield(cfg.BADELEC,'markbad') & cfg.BADELEC.markbad == true
    Marker = union(Marker,cellstr('BAD'));
end

if (isfield(cfg,'STITCH') & ~isempty(cfg.STITCH)) | (isfield(cfg,'STITCHALL') & ~isempty(cfg.STITCHALL))
    Marker = union(Marker,cellstr('Stitch'));
end

if isfield(cfg,'EYEBLINKS') & ~isempty(cfg.EYEBLINKS)
    Marker = union(Marker,cellstr('eyeblink'));
end

if isfield(cfg,'MICROST') & ~isempty(cfg.MICROST)
    if isfield(cfg.MICROST,'doevents_mode') & ...
            strcmp(cfg.MICROST.doevents_mode,'fixed')
        NrClust = cfg.MICROST.doevents;
    else
        NrClust = 20;
    end
    for i = 1:NrClust
        Marker = union(Marker,cellstr(['Micro' num2str(i)]));
    end
end

return