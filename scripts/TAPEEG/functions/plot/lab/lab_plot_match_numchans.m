% Helper file for lab_plot_elec
%
% Written by F. Hatz 2013

function [settings,skipprocessing] = lab_plot_match_numchans(settings,skipdiag)

skipprocessing = 0;
if ~exist('skipdiag','var')
    skipdiag = false;
end
if isempty(settings.DATA)
    return
end
settings = lab_plot_check_numchans(settings);
if isfield(settings,'Match') & settings.Match == true
    return
end
numchans = size(settings.DATA(1).data,2);
if isfield(settings,'LOCS')
    LOCS = settings.LOCS;
    numelectrodes = size(LOCS.x,2);
else
    return
end
if (numelectrodes^2/2 - numelectrodes/2) < numchans
    settings.LOCS = [];
    disp('Mismatch data & loc-information, delete LOCS')
    skipprocessing = 1;
    return
end
numelectrodes = [numelectrodes (numelectrodes^2/2 - numelectrodes/2) (numelectrodes^2/2 + numelectrodes/2)];
if ~isempty(intersect(numelectrodes,numchans))
    settings.Mappings = [];
end
if isfield(settings,'Mappings') & ~isempty(settings.Mappings)
    Mappings = settings.Mappings;
    if (size(Mappings.mappings,2)^2/2 + size(Mappings.mappings,2)/2) >= numchans
        numelectrodes = size(Mappings.mappings,2);
        numelectrodes = [numelectrodes (numelectrodes^2/2 - numelectrodes/2) (numelectrodes^2/2 + numelectrodes/2)];
    else
        disp('Mismatch data & mappings, delete mappings')
        settings.Mappings = [];
        Mappings = [];
    end
else
    Mappings = [];
end
if skipdiag == true & isempty(intersect(numelectrodes,numchans))
    return
end

while skipprocessing == 0 & isempty(intersect(numelectrodes,numchans))
    reduceloc = questdlg('Mismatch loc_file and data! Exclude electrodes or use mappings?', ...
            'Mismatch loc-data','Cancel','Mappings','Exclude','Exclude');
    if strcmp(reduceloc,'Exclude')
        strlist = cellstr(num2str((1:size(LOCS.x,2))'));
        exclude = lab_get_exclude(size(LOCS.x,2));
        if ~isempty(exclude)
            strdefault = setdiff(1:size(LOCS.x,2),exclude);
        else
            strdefault = 1:size(LOCS.x,2);
        end
        selection = listdlg('promptstring','Included channels:','selectionmode','multiple', ...
            'liststring',strlist,'initialvalue',strdefault);
        clearvars strlist strdefault exclude
        if ~isfield(settings,'exclude')
            settings.exclude = setdiff(1:size(LOCS.x,2),selection);
            settings.exclude = settings.exclude(:)';
            settings.numdatachans = size(LOCS.x,2);
        else
            tmp = setdiff(1:size(LOCS.x,2),selection);
            tmp2 = 1:settings.numdatachans;
            tmp2 = setdiff(tmp2,settings.exclude);
            tmp = tmp2(tmp);
            settings.exclude = union(settings.exclude,tmp);
            settings.exclude = settings.exclude(:)';
            clearvars tmp tmp2
        end
        if ~isempty(selection)
            Nchans = [length(selection) (length(selection)^2/2 - length(selection)/2) (length(selection)^2/2 + length(selection)/2)];
            if ~isempty(intersect(Nchans,numchans))
                LOCS = reduce_locs(LOCS,selection);
                numelectrodes = size(LOCS.x,2);
                numelectrodes = [numelectrodes (numelectrodes^2/2 - numelectrodes/2) (numelectrodes^2/2 + numelectrodes/2)]; %#ok<AGROW>
            end
            clearvars Nchans
        else
            skipprocessing = 1;
        end
        clearvars selection
    elseif strcmp(reduceloc,'Mappings')
        if ~isfield(settings,'exclude')
            settings.exclude = [];
            settings.numdatachans = size(LOCS.x,2);
        end
        Mappings = lab_load_mappings([],settings);
        if ~isempty(Mappings) & isstruct(Mappings) & Mappings.mappingsChannels == numelectrodes(1)
            Mchans = length(Mappings.mappings);
            Mchans = [Mchans (Mchans^2/2 - Mchans/2) (Mchans^2/2 + Mchans/2)]; %#ok<AGROW>
        else
            Mchans = [];
        end
        if ~isempty(intersect(Mchans,numchans))
            numelectrodes = length(Mappings.mappings);
            numelectrodes = [numelectrodes (numelectrodes^2/2 - numelectrodes/2) (numelectrodes^2/2 + numelectrodes/2)]; %#ok<AGROW>
        end
    else
        skipprocessing = 1;
    end
end
if skipprocessing == 0
    if numchans == numelectrodes(3)
        [settings.DATA,settings.PLOT] = lab_plot_split_connections(settings.DATA,settings.PLOT,numelectrodes(1));
    end
    settings.LOCS = LOCS;
    if exist('Mappings','var')
        settings.Mappings = Mappings;
        if ~isempty(settings.Mappings) & ~isfield(settings.Mappings,'shortnames')
            shortnames = questdlg('Use short names for display of mappings?', ...
                'Mismatch loc-data','Full','Short','Short');
            if strcmp(shortnames,'Full')
                settings.Mappings.shortnames = false;
            else
                settings.Mappings.shortnames = true;
            end
        end
    end
    settings = lab_plot_check_numchans(settings);
    settings.Match = true;
else
    settings.Match = false;
end

end

function LOCS = reduce_locs(LOCS,selection)
    LOCS.x = LOCS.x(selection);
    LOCS.y = LOCS.y(selection);
    LOCS.z = LOCS.z(selection);
    LOCS.labels = LOCS.labels(selection);
    LOCS.sph_radius = LOCS.sph_radius(selection);
    LOCS.sph_theta = LOCS.sph_theta(selection);
    LOCS.sph_phi = LOCS.sph_phi(selection);
    LOCS.theta = LOCS.theta(selection);
    LOCS.radius = LOCS.radius(selection);
end