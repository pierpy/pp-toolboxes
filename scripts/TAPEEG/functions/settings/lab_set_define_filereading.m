function [cfg,skipprocessing] = lab_set_define_filereading(cfg,header)

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
end

if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var') | ~isfield(cfg,'FILEREADING') | ~isfield(cfg.FILEREADING,'sectionalread')
    cfg.FILEREADING.sectionalread = 0;
end

skipdiag = 1;

Prompt = cell(0,2);
Formats = [];

if ~isfield(header,'EEG_file') | (isfield(header,'EEG_file') & length(header.EEG_file) > 4 & strcmp(header.EEG_file(end-2:end),'mff'))
    Prompt(end+1,:) = {'Select MFF-Segments','SELECTMFF'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_select_mff,'@ALL','@ALL',header};
    skipdiag = 0;
end

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';

if ~isfield(header,'sectionalreading') | header.sectionalreading == true
    Prompt(end+1,:) = {'Sectional Reading (0 = off)','sectionalread'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 80;
    skipdiag = 0;
end

if skipdiag == 0
    [cfg.FILEREADING,Cancelled] = inputsdlg(Prompt,'Reading File',Formats,cfg.FILEREADING);
    if Cancelled == 1 | (cfg.FILEREADING.sectionalread == 0 & isempty(cfg.FILEREADING.SELECTMFF))
        skipprocessing = 1;
        cfg.FILEREADING = [];
    end
else
    cfg.FILEREADING = [];
end

end