function [DATA,PLOT] = lab_plot_split_connections(DATA,PLOT,numchans)

if ~isfield(DATA,'data') | isempty(DATA(1).data)
    return
end
flag = false(size(DATA,1),1);
flag2 = false(size(DATA,1),1);
for i = 1:size(DATA,1)
    if ~isfield(DATA(i,1),'connections') | isempty(DATA(i,1).connections) | DATA(i,1).connections == false
        if isfield(DATA(i,1),'data') & size(DATA(i,1).data,2) == (numchans^2/2 + numchans/2)
            flag(i) = true;
        elseif isfield(DATA(i,1),'data') & size(DATA(i,1).data,2) == (numchans^2/2 - numchans/2)
            flag2(i) = true;
        end
    end
end
if max(flag) == 0
    return
end
flag = repmat(flag,1,size(DATA,2));
flag = flag(:)';
flag2 = repmat(flag2,1,size(DATA,2));
flag2 = flag2(:)';
FlagAll = cat(1,true(1,length(DATA(:))),flag);
FlagAll = FlagAll(:)';
XData = size(DATA,2);
DATA = DATA(:)';
DATA = cat(1,DATA,DATA);
for i = 1:size(DATA,2);
    if flag(i) == true
        [connections,channels] = lab_split_connections(DATA(1,i).data,numchans);
        DATA(1,i).data = connections;
        DATA(2,i).data = channels;
        DATA(1,i).connections = true;
        DATA(2,i).nodes = true;
    elseif flag2(i) == true
        DATA(1,i).connections = true;
    end
end
DATA = DATA(FlagAll);
DATA = reshape(DATA,length(DATA)/XData,XData);
if ~isempty(PLOT)
    XPlot = size(PLOT,2);
    PLOT = PLOT(:)';
    PLOT = cat(1,PLOT,PLOT);
    for i = 1:size(PLOT,2);
        if flag(i) == true
            PLOT(1,i).Name = [PLOT(1,i).Name ' Connections'];
            PLOT(2,i).Name = [PLOT(2,i).Name ' Nodes'];
            PLOT(2,i).AddPlot = true;
        end
    end
    PLOT = PLOT(FlagAll(1:length(PLOT(:))));
    PLOT = reshape(PLOT,length(PLOT)/XPlot,XPlot);
end
        
end