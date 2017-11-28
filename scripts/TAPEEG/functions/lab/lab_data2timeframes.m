function [data,header] = lab_data2timeframes(data,header,TFinclude)
    
header.numtimeframes = length(TFinclude);
if length(TFinclude) ~= size(data,2) & min(TFinclude) >= 1 & max(TFinclude) <= size(data,2)
    tmp = zeros(1,size(data,2));
    tmp(TFinclude) = 1;
    TFinclude = tmp;
    clearvars tmp
end
if length(TFinclude) ~= size(data,2) | min(TFinclude) ~= 0
    return
end
    
if isfield(header,'events') & isfield(header.events,'POS') & ~isempty(header.events.POS)
    header.events.DUR(header.events.DUR<1) = 1;
    events.POS = [];
    events.DUR = [];
    events.OFF = [];
    events.TYP = {};
    for i = 1:length(header.events.POS)
        tmp = zeros(1,size(data,2));
        tmp(header.events.POS(i):header.events.POS(i) + header.events.DUR(i) - 1) = 1;
        tmp = tmp(1,TFinclude==1);
        if max(tmp) == 1
            events.POS = [events.POS int64(find(tmp==1,1,'first'))];
            events.DUR = [events.DUR int64(find(tmp==1,1,'last') - find(tmp==1,1,'first') + 1)];
            if isfield(header.events,'OFF')
                events.OFF = [events.OFF header.events.OFF(i)];
            else
                events.OFF = [events.OFF int64(0)];
            end
            events.TYP = [events.TYP header.events.TYP(i)];
        end
    end
    header.events = events;
end
data = data(:,TFinclude==1);
