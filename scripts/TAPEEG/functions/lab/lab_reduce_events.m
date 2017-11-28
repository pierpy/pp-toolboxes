function events = lab_reduce_events(events,TFinclude)
    
if max(TFinclude) > 1 & min(TFinclude) > 0
    tmp = false(1,max(max(TFinclude),max(events.POS+events.DUR)));
    tmp(TFinclude) = true;
    TFinclude = tmp;
    clearvars tmp
end
    
tmp = find(diff(TFinclude) == -1);
events2.POS = int64(tmp(:)');
events2.DUR = int64(ones(1,length(tmp)));
events2.OFF = int64(zeros(1,length(tmp)));
events2.TYP = repmat({'Cut'},1,length(tmp));
events = lab_mix_markers(events,events2);
clearvars events2 tmp

events.DUR(events.DUR<1) = 1;
events2.POS = [];
events2.DUR = [];
events2.OFF = [];
events2.TYP = {};
Idx = TFinclude == 1;
for i = 1:length(events.POS)
    tmp = zeros(1,length(TFinclude));
    tmp(events.POS(i):events.POS(i) + events.DUR(i) - 1) = 1;
    tmp = tmp(1,Idx);
    if max(tmp) == 1
        events2.POS = [events2.POS int64(find(tmp==1,1,'first'))];
        events2.DUR = [events2.DUR int64(find(tmp==1,1,'last') - find(tmp==1,1,'first') + 1)];
        if isfield(events,'OFF')
            events2.OFF = [events2.OFF events.OFF(i)];
        else
            events2.OFF = [events2.OFF int64(0)];
        end
        events2.TYP = [events2.TYP events.TYP(i)];
    end
end
events = events2;
clearvars events2 tmp

end