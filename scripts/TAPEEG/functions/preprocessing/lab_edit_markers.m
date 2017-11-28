% Edit markers
%
% [data,header,cfg] = lab_edit_markers(data,header,cfg)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
%
% written by F. Hatz 2013

function [data,header,cfg,skipprocessing] = lab_edit_markers(data,header,cfg)
disp('   Edit markers')

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end
if ~exist('header','var') | ~isfield(header,'numchannels')
    header.numchannels = size(data,1);
end

if ~isfield(cfg,'MARK') | ~isfield(cfg.MARK,'edit')
    [cfg,skipprocessing] = lab_set_edit_markers(cfg,header);
    if skipprocessing == 1
        return
    end
end

% Convert old settings
if isfield(cfg.MARK,'markers') & ~isfield(cfg.MARK,'edit')
    cfg.MARK.edit = cfg.MARK.markers;
    cfg.MARK.markers = [];
end

if isfield(cfg,'MARK') & isfield(cfg.MARK,'edit') & ~isempty(cfg.MARK.edit)
    disp('      edit/add markers')
    NewEvents.POS = [];
    NewEvents.DUR = [];
    NewEvents.OFF = [];
    NewEvents.TYP = {};
    if exist('header','var') & isfield(header,'events') & isfield(header.events,'TYP') & ~isempty(header.events.TYP)
        Mkeep = ones(1,length(header.events.POS));
        for i = 1:length(header.events.POS)
            tmp = find(strcmp(cfg.MARK.edit(:,1),header.events.TYP{1,i}));
            if ~isempty(tmp)
                for j = 1:length(tmp)
                    skipPOS = 0;
                    if ~strcmp(cfg.MARK.edit{tmp(j),8},'all')
                        tmp2 = find(strcmp(header.events.TYP,header.events.TYP{1,i}));
                        if strcmp(cfg.MARK.edit{tmp(j),8},'first') & tmp2(1) ~= i
                            skipPOS = 1;
                        elseif strcmp(cfg.MARK.edit{tmp(j),8},'last') & tmp2(end) ~= i
                            skipPOS = 1;
                        end
                        clearvars tmp2
                    end
                    if skipPOS == 0 & cfg.MARK.edit{tmp(j),9} == true
                        Mkeep(1,i) = 0;
                    elseif skipPOS == 0
                        if strcmp(cfg.MARK.edit{tmp(j),3},'seconds') & isfield(header,'samplingrate') & ~isempty(header.samplingrate)
                            shift = cfg.MARK.edit{tmp(j),2} * header.samplingrate;
                            duration = cfg.MARK.edit{tmp(j),4} * header.samplingrate;
                        else
                            shift = cfg.MARK.edit{tmp(j),2};
                            duration = cfg.MARK.edit{tmp(j),4};
                        end
                        if cfg.MARK.edit{tmp(j),6} == true
                            NewEvents.POS = [NewEvents.POS (header.events.POS(1,i) + int64(round(shift)))];
                            if strcmp(cfg.MARK.edit{tmp(j),5},'add')
                                NewEvents.DUR = [NewEvents.DUR (header.events.DUR(1,i) + int64(round(duration)))];
                            else
                                if duration < 0
                                    duration = 0;
                                end
                                NewEvents.DUR = [NewEvents.DUR int64(duration)];
                            end
                            NewEvents.OFF = [NewEvents.OFF int64(0)];
                            if ~isempty(cfg.MARK.edit{tmp(j),7})
                                NewEvents.TYP = [NewEvents.TYP cfg.MARK.edit(tmp(j),7)];
                            else
                                NewEvents.TYP = [NewEvents.TYP cellstr('New Marker')];
                            end
                        else
                            header.events.POS(1,i) = header.events.POS(1,i) + int64(round(shift));
                            if strcmp(cfg.MARK.edit{tmp(j),5},'add')
                                header.events.DUR(1,i) = header.events.DUR(1,i) + int64(round(duration));
                            else
                                if duration < 1
                                    duration = 1;
                                end
                                header.events.DUR(1,i) = int64(round(duration));
                            end
                        end
                    end
                end
            end
        end
        header.events.POS = header.events.POS(1,logical(Mkeep));
        header.events.DUR = header.events.DUR(1,logical(Mkeep));
        header.events.TYP = header.events.TYP(1,logical(Mkeep));
        header.events.OFF = header.events.OFF(1,logical(Mkeep));
    end
    tmp = find(strcmp(cfg.MARK.edit(:,1),'New Marker'));
    for i = 1:length(tmp)
        if strcmp(cfg.MARK.edit{tmp(i),3},'seconds') & isfield(header,'samplingrate') & ~isempty(header.samplingrate)
            shift = cfg.MARK.edit{tmp(i),2} * header.samplingrate;
            duration = cfg.MARK.edit{tmp(i),4} * header.samplingrate;
        else
            shift = cfg.MARK.edit{tmp(i),2};
            duration = cfg.MARK.edit{tmp(i),4};
        end
        if duration < 1
            duration = 1;
        end
        NewEvents.POS = [NewEvents.POS int64(round(shift))];
        NewEvents.DUR = [NewEvents.DUR int64(round(round(duration)))];
        NewEvents.OFF = [NewEvents.OFF int64(0)];
        if ~isempty(cfg.MARK.edit{tmp(i),7})
            NewEvents.TYP = [NewEvents.TYP cfg.MARK.edit(tmp(i),7)];
        else
            NewEvents.TYP = [NewEvents.TYP cellstr('New Marker')];
        end
    end
    if ~isempty(NewEvents.POS)
        header = lab_mix_markers(header,NewEvents);
    end
end

if isfield(cfg,'MARK') & isfield(cfg.MARK,'shuffle') & cfg.MARK.shuffle == true
    disp('      shuffle markers')
    if isfield(cfg.MARK,'markerinclude') & ~isempty(cfg.MARK.markerinclude) & ~strcmp(cfg.MARK.markerinclude{1},'all')
        Markers = cfg.MARK.markerinclude;
    elseif isfield(header,'events') & isfield(header.events,'TYP')
        Markers = unique(header.events.TYP);
    else
        Markers = [];
    end
    if isfield(header,'events')
        header.events = lab_shuffle_markers(header.events,Markers,cfg.MARK.rearrange);
    end
end

if isfield(cfg,'MARK') & isfield(cfg.MARK,'create') & cfg.MARK.create == true
    disp('      create markers')
    header = lab_create_markers(header,cfg.MARK);
end

return