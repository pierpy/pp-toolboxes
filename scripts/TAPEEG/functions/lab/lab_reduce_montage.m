function [montage,include] = lab_reduce_montage(montage,cfg,header,nodiag)

if ~exist('nodiag','var')
    nodiag = false;
end
if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end
if isfield(header,'numdatachannels')
    numchans = header.numdatachannels;
elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans')
    numchans = cfg.EXTRA.numdatachans;
else
    numchans = [];
end

if isnumeric(cfg) & ~isempty(cfg)
    include = cfg;
elseif isfield(header,'includedatachans') & max(header.includedatachans) <= montage(1,1).numchans
    include = header.includedatachans;
elseif ~isempty(numchans) & isfield(cfg,'exclude') & length(setdiff(1:montage(1,1).numchans,cfg.exclude)) == numchans
    include = setdiff(1:montage(1,1).numchans,cfg.exclude);
elseif ~isempty(numchans) & isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude') & ...
        length(setdiff(1:montage(1,1).numchans,cfg.EXTRA.exclude)) == numchans
    include = setdiff(1:montage(1,1).numchans,cfg.EXTRA.exclude);
elseif ~isempty(numchans) & ~isempty(lab_get_exclude(montage(1,1).numchans)) & ...
        length(setdiff(1:montage(1,1).numchans,lab_get_exclude(montage(1,1).numchans))) == numchans
    include = setdiff(1:montage(1,1).numchans,lab_get_exclude(montage(1,1).numchans));
elseif isempty(numchans) & isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    include = 1:montage(1,1).numchans;
    return
elseif isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    include = [];
    montage = [];
    return
elseif nodiag == false
    if isfield(montage(1,1),'includechans') & ~isempty(montage(1,1).includechans)
        strlist = cellstr(num2str(montage(1,1).includechans'));
    else
        strlist = cellstr(num2str((1:montage(1,1).numchans)'));
    end
    if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
        strdefault = cfg.exclude;
    elseif isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'exclude')
        strdefault = cfg.EXTRA.exclude;
    elseif ~isempty(lab_get_exclude(montage(1,1).numchans))
        strdefault = lab_get_exclude(montage(1,1).numchans);
    else
        strdefault = [];
    end
    exclude = listdlg('PromptString','Select excluded channels:','SelectionMode','multiple', ...
        'ListString',strlist,'InitialValue',strdefault,'CancelString','None');
    clearvars strdefault strlist
    include = setdiff((1:montage(1,1).numchans),exclude);
    clearvars exclude
else
    include = 1:montage(1,1).numchans;
    return
end

if isfield(montage(1,1),'includechans') & ~isempty(montage(1,1).includechans)
    include2 = montage(1,1).includechans(include);
else
    include2 = include;
end

keepflag = ones(1,size(montage,2));
for Mnr = 1:size(montage,2)
    maps = zeros(1,montage(1,Mnr).numchans);
    maps(include) = 1:length(include);
    keepflag2 = ones(1,size(montage(1,Mnr).chans,1));
    for i = 1:size(montage(1,Mnr).chans,1)
        if isnumeric(montage(1,Mnr).chans{i,1}) & montage(1,Mnr).chans{i,2} == 0
            montage(1,Mnr).chans{i,1} = maps(montage(1,Mnr).chans{i,1});
            if min(montage(1,Mnr).chans{i,1}) == 0
                keepflag2(1,i) = 0;
            end
        end
        if isnumeric(montage(1,Mnr).chans{i,3}) & ~isnan(montage(1,Mnr).chans{i,3}) & ...
                montage(1,Mnr).chans{i,4} == 0 & montage(1,Mnr).chans{i,3} ~= 0
            montage(1,Mnr).chans{i,3} = maps(montage(1,Mnr).chans{i,3});
            if min(montage(1,Mnr).chans{i,3}) == 0
                keepflag2(1,i) = 0;
            end
        end
    end
    if max(keepflag2) == 1
        montage(1,Mnr).chans = montage(1,Mnr).chans(logical(keepflag2),:);
        montage(1,Mnr).label = montage(1,Mnr).label(logical(keepflag2),:);
    else
        keepflag(1,Mnr) = 0;
    end
    montage(1,Mnr).numchans = length(include);
    montage(1,Mnr).includechans = include2;
end
if max(keepflag) == 1
    montage = montage(1,logical(keepflag));
else
    montage = [];
end