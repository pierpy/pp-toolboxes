function [cfg,skipprocessing] = lab_set_export2edf(cfg,header)

disp ('   EDF-export settings')

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'EDF') | ~isfield(cfg.EDF,'eegsource')
    cfg.EDF.eegsource = 'input';
    cfg.EDF.interpolatebad = true;
    cfg.EDF.filepath = '';
    cfg.EDF.mountpath = '';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Reference','eegsource'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
if ~max(strcmp({'channels','mean','median','laplacian','montage','input'},cfg.EDF.eegsource))
    Formats(end,1).items = {cfg.EDF.eegsource,'channels','mean','median','laplacian','montage','input'};
else
    Formats(end,1).items = {'channels','mean','median','laplacian','montage','input'};
end
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header,'MontageEDF.xls'};
    
Prompt(end+1,:) = {'Montage','montage'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header,'MontageEDF.xls'};

Prompt(end+1,:) = {'Laplacian','LAPL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};

Prompt(end+1,:) = {'Interpolate bad channels','interpolatebad'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 3];

Prompt{end+1,1} = '';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt{end+1,1} = 'Folder to store EDF-file (empty=inputdata-folder)';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'','filepath'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 3];

Prompt{end+1,1} = 'Command to mount storage path (empty = disabled)';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'','mountpath'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).limits = [0 2];
Formats(end,1).size = 250;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Test Command',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@system,[],'mountpath'};

[cfg.EDF,Cancelled] = inputsdlg(Prompt,'EDF-export',Formats,cfg.EDF,3);
pause(0.2);
if isempty(cfg.EDF) | Cancelled == 1
    cfg.EDF = [];
    skipprocessing = 1;
    return
end

end
