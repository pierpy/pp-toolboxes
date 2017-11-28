% Reduce length of EEG/MEG data
%
% [data,header,cfg] = lab_reduce_datalength(data,header,cfg)
%
% data                   = matrix (chans x timeframes)
% header                 = output of lab_read_data
% cfg.REDUCE.firstsample = first sample to include
% cfg.REDUCE.lastsample  = last sample to include
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_reduce_datalength(data,header,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    header = [];
end
if ~isfield(cfg,'REDUCE') | isempty(cfg.REDUCE) | ~isfield(cfg.REDUCE,'firstsample')
    [cfg,skipprocessing] = lab_set_reduce_datalength(cfg,header);
    if skipprocessing == 1 | isempty(cfg.REDUCE.firstsample) | isempty(cfg.REDUCE.lastsample) 
        return
    end
end
newstop = cfg.REDUCE.lastsample;
newstart = cfg.REDUCE.firstsample;
newlength = newstop - newstart + 1;

if newstart == 1 & newstop == size(data,2)
    disp('   Skip reducing length, original length selected')
    return
end

if isfield(cfg,'EDF')
    %--------------------------------------------------------------------------
    % Write first EDF with full datalength (if EDF export is configured)
    %--------------------------------------------------------------------------
    if ~isfield(cfg,'EEG_file') & isfield(header,'EEG_file')
        cfg.EEG_file = header.EEG_file;
        cfg.EEG_filepath = header.EEG_filepath;
    end
    cfg.EDF_file=[cfg.EEG_file(1:end-4) '.edf'];
    cfg.EDF_filepath = cfg.EEG_filepath;
    lab_export2edf(data,header,cfg);
end

%--------------------------------------------------------------------------
% Reduce datalength
%--------------------------------------------------------------------------
if size(data,2) >= newlength
    if size(data,2) < newstop
        newstart = size(data,2) - newlength;
        newstop = size(data,2);
    end
    if newstart == 0
        newstart = 1;
    end
    disp(['Reduce size to TF ' num2str(newstart) '-' num2str(newstop)])
    data = data(:,newstart:newstop);
    if isfield(header,'events') & isfield(header.events,'POS') & ~isempty(header.events.POS)
        header.events.POS =  header.events.POS - int64(newstart-1);
        tmp = find(header.events.POS > 0);
        tmp2 = find(header.events.POS <= size(data,2));
        tmp = intersect(tmp,tmp2);
        header.events.POS =  header.events.POS(1,tmp);
        header.events.TYP =  header.events.TYP(1,tmp);
        header.events.DUR =  header.events.DUR(1,tmp);
        header.events.OFF =  header.events.OFF(1,tmp);
    end
end
header.numtimeframes = (size(data,2));
clearvars tmp tmp2 lag

return