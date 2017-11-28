function [data,header] = lab_microstates2events(data,header,Nr)

disp(['   Convert microstates to events (' num2str(Nr) ' Clusters)'])

if ~isfield(header,'numauxchannels') | header.numauxchannels <= 0 | size(header.channels) < 6
    disp('      Abort: no microstates information found')
    return
end

Mvalues = [];
Cvalues = [];
for i = header.numdatachannels+1:header.numchannels
    if strcmpi(header.channels(i,1:5),'Micro')
        tmp = str2num(header.channels(i,6:end)); %#ok<ST2NM>
        if ~isempty(tmp)
            Mvalues = [Mvalues cat(1,tmp,i)]; %#ok<AGROW>
        end
        clearvars tmp
    elseif strcmpi(header.channels(i,1:4),'Corr')
        tmp = str2num(header.channels(i,5:end)); %#ok<ST2NM>
        if ~isempty(tmp)
            Cvalues = [Cvalues cat(1,tmp,i)]; %#ok<AGROW>
        end
        clearvars tmp
    end
end
if isempty(Mvalues)
    disp('      Abort: no microstates information found')
    return
end

if ~exist('Nr','var') | isempty(Nr)
    Nr = Mvalues(1);
end
tmp = find(Mvalues(1,:)==Nr);
if ~isempty(tmp)
    Mvalues = Mvalues(:,tmp);
else
    disp(['      Abort: microstates information for ' num2str(Nr) ' Clusters not found'])
    return
end
tmp = find(Cvalues(1,:)==Nr);
if ~isempty(tmp)
    Cvalues = Cvalues(2,tmp);
else
    Cvalues = [];
end

Micro = data(Mvalues(2),:);
tmp = find(abs(diff(Micro))>0);
tmp = tmp(:)';
EndMicro = [tmp size(data,2)];
StartMicro = [1 tmp+1];
Idx = find(Micro(StartMicro) ~= 0);
StartMicro = StartMicro(Idx);
EndMicro = EndMicro(Idx);
Ncluster = setdiff(unique(Micro),0);
events.POS = int64(StartMicro(:)');
events.DUR = int64(EndMicro(:)' - StartMicro(:)' + 1);
events.OFF = int64(zeros(1,length(StartMicro)));
events.TYP = cell(1,length(StartMicro));
for i = 1:length(Ncluster)
    Idx = find(Micro(StartMicro) == i);
    events.TYP(1,Idx) = repmat(cellstr(['Micro' num2str(Ncluster(i))]),1,length(Idx));
end
header = lab_mix_markers(header,events);
if ~isempty(Cvalues)
    header.CORR = data(Cvalues,:);
end

end
    