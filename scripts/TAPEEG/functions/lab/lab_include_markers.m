% Helper script to restrict data to included markers
%
% written by F.Hatz 2014

function [data,header,Idx] = lab_include_markers(data,header,markersinclude)

if isfield(header,'events') & ~isempty(header.events)
    events = header.events;
else
    return
end

markers = [];
for i = 1:size(markersinclude,1)
    if strcmp(markersinclude{i,1},'all')
        markers = 1:size(events.POS,2);
        tmp = [];
    else
        tmp = find(strcmp(events.TYP,markersinclude{i,1}));
    end
    markers = union(markers,tmp);
    clearvars tmp
end
Ntf = size(data,2);
TFinclude = zeros(1,Ntf);
for i = 1:length(markers)
    startT = events.POS(1,markers(i));
    endT = events.POS(1,markers(i)) + events.DUR(1,markers(i)) - 1;
    if startT < 1
        startT = 1;
    end
    if endT > Ntf;
        endT = Ntf;
    end
    if startT < endT
        TFinclude(1,startT:endT) = 1;
    else
        disp(['      skip marker ' events.TYP{1,markers(i)} ', not valid or duration <= 1'])
    end
end
if max(TFinclude) == 0
    data = [];
elseif min(TFinclude) == 0
    Idx = find(TFinclude == 1);
    header.events = lab_reduce_events(events,TFinclude);
    data = data(:,Idx);
    header.numtimeframes = size(data,2);
end
clearvars markers i

end