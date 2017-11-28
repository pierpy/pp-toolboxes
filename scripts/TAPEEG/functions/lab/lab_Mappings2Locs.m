function LOCSout = lab_Mappings2Locs(Mappings,LOCS)

LOCSout = [];
if isempty(Mappings)  | ~isfield(Mappings,'mappingsChannels') | isempty(LOCS) |  ...
        ~isfield(LOCS,'x') | Mappings.mappingsChannels ~= size(LOCS.x,2)
    return
end

LOCSout.aux = LOCS.aux;
LOCSout.auxlabels = LOCS.auxlabels;

for i = 1:size(Mappings.mappings,2)
    LOCSout.x(1,i) = mean(LOCS.x(1,Mappings.mappings{i}));
    LOCSout.y(1,i) = mean(LOCS.y(1,Mappings.mappings{i}));
    LOCSout.z(1,i) = mean(LOCS.z(1,Mappings.mappings{i}));
end

if isfield(Mappings,'shortnames') & Mappings.shortnames == true
    LOCSout.labels = Mappings.mappingstitleS(:)';
else
    LOCSout.labels = Mappings.mappingstitle(:)';
end

LOCSout = lab_locs2sph(LOCSout);

return