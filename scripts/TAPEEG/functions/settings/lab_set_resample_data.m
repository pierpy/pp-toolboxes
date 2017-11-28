function [cfg,skipprocessing] = lab_set_resample_data(cfg,header)
    
skipprocessing = 0;

disp ('   Ask for new samplingrate')

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if  ~exist('cfg','var') | ~isfield(cfg,'RESAMPLE') | ~isfield(cfg.RESAMPLE,'new_samplingrate')
    if exist('header','var') & isfield(header,'samplingrate')
        cfg.RESAMPLE.new_samplingrate = header.samplingrate;
    else
        cfg.RESAMPLE.new_samplingrate = [];
    end
end

numel = 1;
Prompt(numel,:) = {'New samplingrate','new_samplingrate'};
Formats(numel,1).type = 'edit';
Formats(numel,1).format = 'float';
Formats(numel,1).limits = [0 inf]; % 2-digits (positive #)
Formats(numel,1).size = 45;

[cfg.RESAMPLE,Cancelled] = inputsdlg(Prompt,'Re-sampling',Formats,cfg.RESAMPLE);
if ~isfield(cfg.RESAMPLE,'new_samplingrate') | isempty(cfg.RESAMPLE.new_samplingrate) | Cancelled == 1
    cfg.RESAMPLE = [];
    skipprocessing = 1;
    return
end