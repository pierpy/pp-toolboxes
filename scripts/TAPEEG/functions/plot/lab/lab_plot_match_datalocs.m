% Helper file for lab_plot_elec
%
% Written by F. Hatz 2013

function [settings,skipprocessing] = lab_plot_match_datalocs(settings,numchans,nomapping)

skipprocessing = 0;

if ~exist('nomapping','var')
    nomapping = false;
end
if ~exist('numchans','var')
    answer = inputdlg({'number of channels in data'},'Number channels',[1 50]);
    numchans = str2num(answer{1,1}); %#ok<ST2NM>
    clearvars answer
end
if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'doconnections')
    settings.doconnections = false;
end

if isstruct(settings) & isfield(settings,'LOCS')
    LOCS = settings.LOCS;
elseif isstruct(settings) & isfield(settings,'x')
    LOCS = settings;
    settings = [];
end
if exist('LOCS','var') & ischar(LOCS) & exist(LOCS,'file')
    LOCS = lab_read_locs(LOCS);
end
if ~exist('LOCS','var') | ~isstruct(LOCS)
    LOCS = lab_read_locs;
end
numelectrodes = size(LOCS.x,2);
if nomapping == false
    if isstruct(settings) & isfield(settings,'mappings') & ~isempty(settings.mappings)
        mappings = settings.mappings;
        numelectrodes = size(mappings.mappings,2);
    end
    if exist('mappings','var') & ischar(mappings) & exist(mappings,'file')
        mappings = lab_read_mappings(mappings);
        numelectrodes = size(mappings.mappings,2);
    end
end
if isempty(numchans) | numchans == 1
    settings.LOCS = LOCS;
    if exist('mappings','var')
        settings.mappings = mappings;
    end
    return
end

if numchans == numelectrodes^2/2 + numelectrodes/2
    settings.doconnections = true;
    numelectrodes = numchans;
elseif numchans == numelectrodes^2/2 - numelectrodes/2
    settings.doconnections = true;
    numelectrodes = numchans;
end

if numelectrodes < numchans
    disp('Abort: Mismatch data & loc-information')
    skipprocessing = 1;
    return
end

while skipprocessing == 0 & numelectrodes > numchans
    if nomapping == false
        reduceloc = questdlg('Mismatch loc_file and data! Exclude electrodes or use mappings?', ...
            'Mismatch loc-data','Cancel','Mappings','Exclude','Exclude');
    else
        reduceloc = questdlg('First match data to electrode_file by excluding electrodes!','Mismatch loc-data','Cancel','Exclude','Exclude');
    end
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
        if ~isempty(selection) & numchans == length(selection)
            LOCS = reduce_locs(LOCS,selection);
            numelectrodes = size(LOCS.x,2);
        elseif ~isempty(selection) & numchans == length(selection)^2/2 + length(selection)/2
            LOCS = reduce_locs(LOCS,selection);
            numelectrodes = length(selection)^2/2 + length(selection)/2;
            settings.doconnections = true;
        elseif ~isempty(selection) & numchans == length(selection)^2/2 - length(selection)/2
            LOCS = reduce_locs(LOCS,selection);
            numelectrodes = length(selection)^2/2 - length(selection)/2;
            settings.doconnections = true;
        elseif isempty(selection)
            skipprocessing = 1;
        end
        clearvars selection
    elseif strcmp(reduceloc,'Mappings')
        if ~isfield(settings,'exclude')
            settings.exclude = [];
            settings.numdatachans = size(LOCS.x,2);
        end
        mappings = lab_load_mappings([],settings);
        if ~isempty(mappings) & isstruct(mappings) & mappings.mappingsChannels == numelectrodes
            Mchans = length(mappings.mappings);
        else
            Mchans = [];
        end
        if numchans == Mchans
            numelectrodes = Mchans;
        elseif numchans == Mchans^2/2 + Mchans/2
            numelectrodes = Mchans^2/2 + Mchans/2;
            settings.doconnections = true;
        elseif numchans == Mchans^2/2 - Mchans/2
            numelectrodes = Mchans^2/2 - Mchans/2;
            settings.doconnections = true;
        end
    else
        skipprocessing = 1;
    end
end
if skipprocessing == 0
    settings.LOCS = LOCS;
    if exist('mappings','var')
        settings.mappings = mappings;
        if ~isempty(settings.mappings) & ~isfield(settings,'shortnames')
            shortnames = questdlg('Use short names for display of mappings?', ...
                'Mismatch loc-data','Full','Short','Short');
            if strcmp(shortnames,'Full')
                settings.shortnames = false;
            else
                settings.shortnames = true;
            end
        end
    end
end
pause(0.2);

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