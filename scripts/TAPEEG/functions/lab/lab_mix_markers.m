function header = lab_mix_markers(header,headertmp)

if isfield(header,'POS')
    header.events = header;
    rmheader = 1;
else
    rmheader = 0;
end
if isfield(headertmp,'POS')
    headertmp.events = headertmp;
end

if isfield(headertmp,'events') & ~isempty(headertmp.events) & ...
        isfield(headertmp.events,'POS') & ~isempty(headertmp.events.POS)
    if ~isfield(header,'events') | ~isfield(header.events,'POS')
        header.events.POS = [];
        header.events.DUR = [];
        header.events.TYP = [];
        header.events.OFF = [];
    end
    header.events.POS = [header.events.POS headertmp.events.POS];
    header.events.DUR = [header.events.DUR headertmp.events.DUR];
    header.events.TYP = [header.events.TYP headertmp.events.TYP];
    header.events.OFF = [header.events.OFF headertmp.events.OFF];
    tmp = 1:size(header.events.POS,2);
    tmp(2,:) = header.events.POS;
    tmp = sortrows(tmp',2)';
    tmp = tmp(1,:);
    header.events.POS = header.events.POS(tmp);
    header.events.DUR = header.events.DUR(tmp);
    header.events.TYP = header.events.TYP(tmp);
    header.events.OFF = header.events.OFF(tmp);
end

if isfield(header,'events')& isfield(header.events,'POS') & ~isempty(header.events.POS) 
    % delete doubles
    [~,P1] = unique(header.events.POS,'first');
    [~,P2] = unique(header.events.POS,'last');
    Idx = true(1,length(header.events.POS));
    for i = 1:length(P1)
        if P1(i) ~= P2(i)
            for j = P1(i)+1:P2(i)
                if header.events.DUR(j) == header.events.DUR(j-1) & strcmp(header.events.TYP{j},header.events.TYP{j-1})
                    Idx(j) = false;
                end
            end
        end
    end
    header.events.POS = header.events.POS(1,Idx);
    header.events.DUR = header.events.DUR(1,Idx);
    header.events.OFF = header.events.OFF(1,Idx);
    header.events.TYP = header.events.TYP(1,Idx);
end

if rmheader == 1
    header = header.events;
end

return