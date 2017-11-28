function settings = lab_plot_check_numchans(settings)

settings.Match = false;
if isempty(settings.DATA)
    return
end
if isfield(settings,'LOCS') & ~isempty(settings.LOCS)
    numelectrodes = size(settings.LOCS.x,2);
elseif isfield(settings,'Brain') & isfield(settings.Brain,'labels')
    numelectrodes = length(settings.Brain.labels);
else
    if ~isempty(settings.PLOT)
        for i = 1:size(settings.PLOT,1)
            for j = 1:size(settings.PLOT,2)
                settings.PLOT(i,j).Valid = false;
            end
        end
    end
    return
end
if isfield(settings,'Mappings') & ~isempty(settings.Mappings)
    numelectrodes = size(settings.Mappings.mappings,2);
end
if isempty(settings.PLOT)
    for i = 1:size(settings.DATA,1)
        for j = 1:size(settings.DATA,2)
            numchans = size(settings.DATA(i,j).data,2);
            if numchans == numelectrodes
                settings.Match = true;
            elseif numchans == (numelectrodes^2/2 - numelectrodes/2)
                settings.Match = true;
            end
        end
    end
elseif size(settings.PLOT,2) == 1
    for i = 1:size(settings.DATA,1)
        numchans = size(settings.DATA(i,1).data,2);
        if numchans == numelectrodes
            if ~isfield(settings.PLOT,'Mode') | length(settings.PLOT) < i | ~any(strcmp({'Nodes','Surface','Volume'},settings.PLOT(i,1).Mode))
                settings.PLOT(i,1).Mode = 'Nodes';
            end
            settings.PLOT(i,1).Valid = true;
            settings.Match = true;
        elseif numchans == (numelectrodes^2/2 - numelectrodes/2)
            settings.PLOT(i,1).Mode = 'Connections';
            settings.PLOT(i,1).Valid = true;
            settings.Match = true;
        else
            settings.PLOT(i,1).Valid = false;
        end
    end
else
    for i = 1:size(settings.DATA,1)
        for j = 1:size(settings.DATA,2)
            numchans = size(settings.DATA(i,j).data,2);
            if numchans == numelectrodes
                if ~isfield(settings.PLOT,'Mode') | ~any(strcmp({'Nodes','Surface','Volume'},settings.PLOT(i,j).Mode))
                    settings.PLOT(i,j).Mode = 'Nodes';
                end
                settings.PLOT(i,1).Valid = true;
                settings.Match = true;
            elseif numchans == (numelectrodes^2/2 - numelectrodes/2)
                settings.PLOT(i,j).Mode = 'Connections';
                settings.PLOT(i,1).Valid = true;
                settings.Match = true;
            else
                settings.PLOT(i,1).Valid = false;
            end
        end
    end
end

end