function elec = lab_locs2elec(LOCS)

if isfield(LOCS,'x')
    elec.chanpos(:,1) = LOCS.x';
    elec.chanpos(:,2) = LOCS.y';
    elec.chanpos(:,3) = LOCS.z';
    elec.elecpos = elec.chanpos;
    elec.label = LOCS.labels';
    elec.unit = 'mm';
else
    elec = [];
end

return