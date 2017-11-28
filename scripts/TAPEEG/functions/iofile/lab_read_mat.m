% Function to read Netstation export (.mat)
%
% [data,header,cfg] = lab_read_mat(filename,cfg)
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_mat(Filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end

MAT = load(Filename);
if isfield(MAT,'Result') | (isfield(MAT,'result') & isfield(MAT,'patient'))
    [data,header,cfg] = lab_read_matrix(Filename,cfg);
    return
end

Vars = fieldnames(MAT);
if isfield(MAT,'Impedances_0')
    nChans = size(MAT.Impedances_0,1);
else
    nChans = [];
    for i = 1:size(Vars,1);
        if isempty(nChans) | size(MAT.(Vars{i,1}),1) > nChans
            nChans = size(MAT.(Vars{i,1}),1);
        end
    end
end
namedata = [];
nameextdata = [];
numberextdata = [];
numtimeframes = [];
events = [];
for i = 1:size(Vars,1);
    if isa(MAT.(Vars{i,1}),'double') & size(MAT.(Vars{i,1}),2) > 1 & ~isempty(strfind(Vars{i,1},'_'))
        if size(MAT.(Vars{i,1}),1) == nChans
            numtimeframes(1,end+1) = size(MAT.(Vars{i,1}),2); %#ok<AGROW>
            namedata{1,end+1} = Vars{i,1}; %#ok<AGROW>
            nameextdata{1,end+1} = ''; %#ok<AGROW>
            numberextdata(1,end+1) = 0; %#ok<AGROW>
        else
            tmp = [];
            for j = 1:length(namedata)
                if ~isempty(strfind(Vars{i,1},namedata{1,j}))
                    tmp = j;
                end
            end
            if ~isempty(tmp)
                nameextdata{1,tmp(1)} = Vars{i,1}; %#ok<AGROW>
                numberextdata(1,tmp(1)) = size(MAT.(Vars{i,1}),1); %#ok<AGROW>
            end
        end
    end
    if size(MAT.(Vars{i,1}),1) == 4 & isa(MAT.(Vars{i,1}),'cell')
        events = [events Vars(i,1)]; %#ok<AGROW>
    end
end
clearvars i j tmp

if ~isempty(numtimeframes)
    % Write header
    header.channels = num2str((1:nChans)');
    numberextdata = min(numberextdata);
    if numberextdata > 0
        for i = 1:numberextdata
            tmp = ['e' num2str(i)];
            header.channels(nChans+i,1:length(tmp)) = tmp;
        end
    else
        numberextdata = 0;
    end
    tmp = strfind(namedata{1,1},'_');
    header.version=0;
    header.numchannels=nChans;
    header.numauxchannels=numberextdata;
    header.numtimeframes=sum(numtimeframes);
    if exist('samplingRate','var')
        header.samplingrate=samplingRate;
    else
        header.samplingrate = [];
    end
    header.subject.name = namedata{1,1}(1:(tmp(1)-1));
    string = namedata{1,1}((tmp(end-1)+1):(tmp(end)-1));
    if length(string) == 8
        header.year = string(1:4);
        header.month = string(5:6);
        header.day = string(7:8);
    elseif length(string) == 4
        header.year = 1900;
        header.month = string(1:2);
        header.day = string(3:4);
    end
    string = namedata{1,1}((tmp(end)+1):end);
    header.hour = string(1:2);
    header.minute = string(3:4);
    header.second=0;
    header.millisecond=0;
    data = [];
    for i = 1:size(namedata,2)
        if ~isempty(nameextdata{1,i})
            datatmp = cat(1,(MAT.(namedata{1,i})),(MAT.(nameextdata{1,i})));
        else
            datatmp = MAT.(namedata{1,i});
        end
        data = [data datatmp(1:(header.numchannels+header.numauxchannels),:)]; %#ok<AGROW>
    end
    clearvars datatmp tmp
    
    % Calculate events
    offset = 0;
    if size(numtimeframes,2) > 1
        for i = 2:size(numtimeframes,2)
            offset(1,i) = sum(numtimeframes(1:(i-1))); %#ok<AGROW>
        end
    end
    event.TYP = [];
    event.POS = [];
    event.DUR = [];
    event.OFF = [];
    for i = 1:size(events,2)
        tmp = MAT.(events{1,i});
        for j = 1:size(tmp,2)
            event.TYP = [event.TYP tmp(1,j)];
            event.POS = [event.POS int64(tmp{2,j} + offset(1,tmp{3,j}))];
            event.DUR = [event.DUR int64(tmp{4,j} - tmp{2,j} + 1)];
            event.OFF = [event.OFF int64(0)];
        end
    end
    if ~isempty(event.POS)
        tmp = cell(4,size(event.POS,2));
        tmp(1,:) = event.TYP;
        tmp(2,:) = num2cell(event.POS);
        tmp(3,:) = num2cell(event.DUR);
        tmp(4,:) = num2cell(event.OFF);
        tmp = sortrows(tmp',2)';
        header.events.TYP = tmp(1,:);
        header.events.POS = cell2mat(tmp(2,:));
        header.events.DUR = cell2mat(tmp(3,:));
        header.events.OFF = cell2mat(tmp(4,:));
    end
    clearvars event
else
    data = [];
    header = [];
end

return