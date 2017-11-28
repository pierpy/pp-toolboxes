function header = lab_reallocate_markers(header,timeline,markerexclude)

if ~exist('timeline','var')
    timeline = 1:header.numtimeframes;
end
if ~isfield(header,'events') | ~isfield(header.events,'POS') | isempty(header.events.POS)
    return
end

tmp = timeline<=header.numtimeframes;
timeline = timeline(tmp);

if exist('markerexclude','var') & ~isempty(markerexclude)
    Sel = [];
    for i = 1:length(markerexclude)
        Sel = union(tmp,find(ismember(header.events.TYP,markerexclude{i})==1));
    end
    Sel = setdiff(1:length(header.events.POS),Sel);
else
    Sel = 1:length(header.events.POS);
end

events.POS = [];
events.OFF = [];
events.DUR = [];
events.TYP = {};
for i = Sel
    if header.events.DUR(i) > 1
        tmp = zeros(1,header.numtimeframes);
        tmp(header.events.POS(i):header.events.POS(i)+header.events.DUR(i)-1) = 1;
        tmp = tmp(timeline);
        startE = find(tmp==1,1,'first');
        stopE = find(tmp==1,1,'last');
        if ~isempty(startE)
            events.POS(1,end+1) = int64(startE);
            events.DUR(1,end+1) = int64(stopE) - int64(startE) + 1;
            events.DUR(1,end+1) = int64(0);
            events.TYP{1,end+1} = header.events.TYP{i};
        end
    else
        tmp = find(timeline == header.events.POS(i));
        if ~isempty(tmp)
            events.POS(1,end+1) = int64(tmp(1));
            events.DUR(1,end+1) = int64(1);
            events.DUR(1,end+1) = int64(0);
            events.TYP{1,end+1} = header.events.TYP{i};
        end
    end
end
if ~isempty(events.POS)
    events.POS = int64(events.POS);
    events.OFF = int64(events.OFF);
    events.DUR = int64(events.DUR);
    header.events = events;
else
    header.events = [];
end
header.numtimeframes = length(timeline);