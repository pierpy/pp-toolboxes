function [cfg,skipprocessing] = lab_set_reduce_datalength(cfg,header)

skipprocessing = 0;

if ~exist('header','var')
    header =[];
end
if ~exist('cfg','var') | ~isfield(cfg,'REDUCE') | ~isfield(cfg.REDUCE,'firstsample')
    cfg.REDUCE.firstsample = 1;
end
if isfield(header,'numtimeframes') & ~isfield(cfg.REDUCE,'lastsample')
    cfg.REDUCE.lastsample = header.numtimeframes;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'First sample','firstsample'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999];
Formats(end,1).size = 50;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Last sample','lastsample'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999];
Formats(end,1).size = 50;

[cfg.REDUCE,Cancelled] = inputsdlg(Prompt,'Reduce Length',Formats,cfg.REDUCE);
if Cancelled == 1
    cfg.REDUCE = [];
    skipprocessing = 1;
end

return