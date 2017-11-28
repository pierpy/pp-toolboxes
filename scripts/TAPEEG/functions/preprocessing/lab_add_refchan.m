function [header,data,cfg] = lab_add_refchan(header,data,cfg,novrb)
    
if ~exist('novrb','var')
    novrb = false;
end
if novrb == false
    disp('Add reference channel (empty channel)')
end

if isfield(header,'ref_chan') & isnumeric(header.ref_chan) & ~isempty(header.ref_chan)
    disp('   Abort: Ref-Channel already existing')
    return
end

if header.numdatachannels < header.numchannels
    auxtmp = data(header.numdatachannels+1:end,:);
    auxchannels = cellstr(header.channels(header.numdatachannels+1:end,:));
    if isfield(header,'eog_ch') & header.eog_ch > header.numdatachannels
        header.eog_ch = header.eog_ch + 1;
    end
    if isfield(header,'ecg_ch') & header.ecg_ch > header.numdatachannels
        header.ecg_ch = header.ecg_ch + 1;
    end
else
    auxtmp = [];
    auxchannels = {};
end
channels = cellstr(header.channels(1:header.numdatachannels,:));

if ~isempty(data)
    data = data(1:header.numdatachannels,:);
    data(end+1,:) = 0;
    if ~isempty(auxtmp)
        data = cat(1,data,auxtmp);
    end
    header.numchannels = size(data,1);
else
    header.numchannels = header.numchannels + 1;
end
header.numdatachannels = header.numdatachannels + 1;
header.ref_chan = header.numdatachannels;


if isfield(cfg,'ADDREF') & isfield(cfg.ADDREF,'name') & ~isempty(cfg.ADDREF.name)
    channels{end+1,1} = cfg.ADDREF.name;
else
    channels{end+1,1} = 'REF';
end
header.channels = char(cat(1,channels,auxchannels));

if isfield(header,'locs') & isfield(header.locs,'x') & length(header.locs.x) == header.numdatachannels-1
    if ~strcmp(channels{end},'REF')& isfield(cfg,'LOCS') & ...
            isfield(cfg.LOCS,'matchlocs') & cfg.LOCS.matchlocs == true
        Idx = find(strcmp(cfg.LOCS.locs.labels,channels{end}));
        if ~isempty(Idx)
            header.locs.x(1,end+1) = cfg.LOCS.locs.x(Idx(1));
            header.locs.y(1,end+1) = cfg.LOCS.locs.y(Idx(1));
            header.locs.z(1,end+1) = cfg.LOCS.locs.z(Idx(1));
            header.locs.labels{1,end+1} = channels{end};
            header.locs = lab_locs2sph(header.locs);
        else
            header.locs.x(1,end+1) = 0;
            header.locs.y(1,end+1) = 0;
            header.locs.z(1,end+1) = 0;
            header.locs.labels{1,end+1} = channels{end};
            header.locs = lab_locs2sph(header.locs);
        end
    else
        header.locs.x(1,end+1) = 0;
        header.locs.y(1,end+1) = 0;
        header.locs.z(1,end+1) = 0;
        header.locs.labels{1,end+1} = channels{end};
        header.locs = lab_locs2sph(header.locs);
    end
end

end