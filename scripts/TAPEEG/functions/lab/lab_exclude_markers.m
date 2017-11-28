% Helper script to exclude markers in data file
%
% written by F.Hatz 2014

function [data,header,Idx] = lab_exclude_markers(data,header,markerexclude)

if isfield(header,'events') & ~isempty(header.events)
    events = header.events;
else
    Idx = 1:size(data,2);
    return
end

markers = [];
for i = 1:size(markerexclude,1)
    if strcmp(markerexclude{i,1},'all')
        markers = 1:size(events.POS,2);
        tmp = [];
    else
        tmp = find(strcmp(events.TYP,markerexclude{i,1}));
    end
    markers = union(markers,tmp);
    clearvars tmp
end
Ntf = size(data,2);
TFinclude = ones(1,Ntf);
for i = 1:length(markers)
    startT = events.POS(1,markers(i));
    endT = events.POS(1,markers(i)) + events.DUR(1,markers(i)) - 1;
    if startT < 1
        startT = 1;
    end
    if endT > Ntf
        endT = Ntf;
    end
    if startT < endT
        TFinclude(1,startT:endT) = 0;
    else
        disp(['      skip marker ' events.TYP{1,markers(i)} ', not valid or duration <= 1'])
    end
end
if max(TFinclude) == 0
    data = [];
    Idx = [];
elseif min(TFinclude) == 0
    Idx = find(TFinclude == 1);
    header.events = lab_reduce_events(events,TFinclude);
    data = data(:,Idx);
    header.numtimeframes = size(data,2);
else
    Idx = 1:size(data,2);
end
clearvars markers i

end