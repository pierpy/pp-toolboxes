function [cfg,skipprocessing] = lab_set_references(cfg,header)
disp ('   Ask for references-settings')
skipprocessing = 0;

if ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'REF') | isempty(cfg.REF)
    cfg.REF.eegsource = {'mean'};
    cfg.REF.LAPL.lap_maxdistance=4;
    cfg.REF.LAPL.lap_weightmaxdistance=5;
    cfg.REF.LAPL.lap_excluderef = true;
    cfg.REF.montage = [];
end

for i = 1:length(cfg.REF.eegsource)
    refitems = {'mean';'median';'laplacian';'channels';'montage'};
    if ~max(strcmp({'mean','median','laplacian','channels','montage'},cfg.REF.eegsource{i}))
        refitems = cat(1,refitems,cfg.REF.eegsource{i});
    end
end

Prompt = {};
Formats = {};

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Reference','eegsource'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = refitems;
Formats(end,1).size = [80 90];
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};
Formats(end,1).span = [2 1];

Prompt(end+1,:) = {'Montage','montage'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header};

Prompt(end+1,:) = {'Laplacian','LAPL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};

[cfg.REF,Cancelled] = inputsdlg(Prompt,'Reference data',Formats,cfg.REF);
if isempty(cfg.REF) | Cancelled == 1
    cfg.REF = [];
    skipprocessing = 1;
    return
end