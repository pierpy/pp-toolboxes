function [DATA2,PLOT2] = lab_plot_create_title(DATA2,PLOT2)

Modes = fieldnames(DATA2);
for i = 1:length(DATA2)
    Measure = {};
    Subject = {};
    for j = 1:length(Modes)
        if ~isempty(DATA2(i).(Modes{j})) & ~any(strcmp(Measure,DATA2(i).(Modes{j}).measure))
            Measure{end+1} = DATA2(i).(Modes{j}).measure; %#ok<AGROW>
        end
        if ~isempty(DATA2(i).(Modes{j})) & ~any(strcmp(Subject,DATA2(i).(Modes{j}).subject))
            Subject{end+1} = DATA2(i).(Modes{j}).subject; %#ok<AGROW>
        end
    end
    if ~isempty(Subject) & ~isempty(Subject{1}) & ~isempty(Measure) & ~isempty(Measure{1})
        PLOT2(i).Name = [sprintf(Measure{:},'%s_') '_' sprintf(Subject{:},'%s_')];
    elseif ~isempty(Subject) & ~isempty(Subject{1})
        PLOT2(i).Name = sprintf(Subject{:},'%s_');
    elseif ~isempty(Measure) & ~isempty(Measure{1})
        PLOT2(i).Name = sprintf(Measure{:},'%s_');
    else
        PLOT2(i).Name = 'PLOT';
    end
end

end