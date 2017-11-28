function header = lab_create_markers(header,settings)
    
if ~exist('settings','var') | isempty(settings);
    settings = [];
    settings = lab_set_create_markers(settings);
end
markers = settings.markers;
if strcmp(markers{1},'all')
    if ~isfield(header,'events') & isfield(header.events,'TYP')
        markers = unique(header.events.TYP);
    else
        return
    end
end

if isfield(header,'events') & ~isempty(header.events)
    events = header.events;
    tmp = unique(events.TYP);
    EventsLength = cell(2,0);
    for i = 1:length(tmp);
        Idx = strcmp(events.TYP,tmp{i});
        EventsLength{1,i} = tmp{i};
        EventsLength{2,i} = mean(events.DUR(1,Idx)+1);
    end
    clearvars tmp i Idx
    Idx = [];
    for i = 1:length(markers)
        tmp = find(strcmp(events.TYP,markers{i}));
        Idx = [Idx tmp(:)']; %#ok<AGROW>
    end
    Idx = setdiff(1:length(events.TYP),Idx);
    events.POS = events.POS(1,Idx);
    events.DUR = events.DUR(1,Idx);
    events.OFF = events.OFF(1,Idx);
    events.TYP = events.TYP(1,Idx);
    clearvars Idx
else
    events.POS = [];
    events.DUR = [];
    events.OFF = [];
    events.TYP = {};
    EventsLength = {' ';[]};
end

events2.POS = int64(1);
tmp = find(strcmp(EventsLength(1,:),markers{1}), 1);
if ~isempty(tmp)
    events2.DUR = int64(ceil(EventsLength{2,tmp}));
elseif isfield(settings.duration) & ~isempty(settings.duration)
    events2.DUR = int64(ceil(settings.duration));
else
    events2.DUR = int64(1);
end
events2.OFF = int64(0);
events2.TYP = markers(1);
for i = 2:length(markers)
    events2.POS(end+1) = events2.POS(1,i-1) + events2.DUR(1,i-1) + 1;
    tmp = find(strcmp(EventsLength(1,:),markers{i}), 1);
    if ~isempty(tmp)
        events2.DUR(1,end+1) = int64(ceil(EventsLength{2,tmp}));
    elseif isfield(settings.duration) & ~isempty(settings.duration)
        events2.DUR(1,end+1) = int64(ceil(settings.duration));
    else
        events2.DUR(1,end+1) = int64(1);
    end
    events2.OFF(1,end+1) = int64(0);
    events2.TYP(1,end+1) = markers(i);
end
L = events2.POS(end) + events2.DUR(end);
NumRep = ceil(header.numtimeframes / L);
Factor = repmat(L*(0:NumRep-1),length(events2.POS),1);
Factor = Factor(:)';
events2.POS = repmat(events2.POS,1,NumRep) + Factor;
events2.DUR = repmat(events2.DUR,1,NumRep);
events2.OFF = repmat(events2.OFF,1,NumRep);
events2.TYP = repmat(events2.TYP,1,NumRep);
Idx = (events2.POS + events2.DUR) < header.numtimeframes;
events2.POS = events2.POS(1,Idx);
events2.DUR = events2.DUR(1,Idx);
events2.OFF = events2.OFF(1,Idx);
events2.TYP = events2.TYP(1,Idx);
header.events = lab_mix_markers(events,events2);

return