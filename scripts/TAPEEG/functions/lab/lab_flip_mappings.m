% If a mappings is only defined for left or right hemisphere, this function
% flips the mapping to the other side
%
% Mappings = lab_flip_mappings(Mappings,Locs)
%
% written by F. Hatz 2014

function Mappings = lab_flip_mappings(Mappings,Locs)

if ~exist('Mappings','var')
    [Mappings,Mappings_file] = lab_load_mappings;
elseif ischar(Mappings)
    Mappings_file = Mappings;
    Mappings = lab_read_data(Mappings);
end
if ~exist('Locs','var')
    Locs = lab_read_locs;
elseif ischar(Locs)
    Locs = lab_read_locs(Locs);
end
if isempty(Mappings) | ~isstruct(Mappings) | ~isfield(Mappings,'mappings') | ...
        isempty(Locs) | ~isfield(Locs,'x')
    return
end

if Mappings.mappingsChannels ~= length(Locs.x)
    return
end
Locs = [Locs.x' Locs.y' Locs.z'];

Nmaps = length(Mappings.mappings);
for i = 1:Nmaps
    selection = Mappings.mappings{i};
    selection2 = [];
    for j = 1:length(selection)
        loctmp = Locs(selection(j),:);
        loctmp(1) = -loctmp(1);
        distance = lab_distance(Locs,loctmp);
        selection2 = [selection2 find(distance == min(distance),1)];
    end
    Mappings.mappings{1,end+1} = selection2;
    if strcmp(Mappings.mappingstitle{i,1}(1),'L')
        Mappings.mappingstitle{end+1,1} = ['R' Mappings.mappingstitle{i,1}(2:end)];
    elseif strcmp(Mappings.mappingstitle{i,1}(1),'R')
        Mappings.mappingstitle{end+1,1} = ['L' Mappings.mappingstitle{i,1}(2:end)];
    elseif strcmp(Mappings.mappingstitle{i,1}(end),'R')
        Mappings.mappingstitle{end+1,1} = [Mappings.mappingstitle{i,1}(1:end-1) 'L'];
    elseif strcmp(Mappings.mappingstitle{i,1}(end),'L')
        Mappings.mappingstitle{end+1,1} = [Mappings.mappingstitle{i,1}(1:end-1) 'R'];
    else
        Mappings.mappingstitle{end+1,1} = [Mappings.mappingstitle{i,1} '_flip'];
    end
end

if exist('Mappings_file','var')
    [~,Mappings_filepath,~,Mappings_fileS] = lab_filename(Mappings_file)
    lab_write_mapping(fullfile(Mappings_filepath,[Mappings_fileS '_flip.xls']),Mappings);
end