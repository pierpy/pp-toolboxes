function [cfg,skipprocessing] = lab_set_addref(cfg)

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end

if ~isfield(cfg,'ADDREF') | ~isfield(cfg.ADDREF,'name') | isempty(cfg.ADDREF.name)
    cfg.ADDREF.name = 'REF';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Ref-Channel','name'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 50;
        
[cfg.ADDREF,Cancelled] = inputsdlg(Prompt,'Extra channels',Formats,cfg.ADDREF);
pause(0.1);
if Cancelled == 1
    cfg.ADDREF = [];
    skipprocessing = 1;
    return
end

end