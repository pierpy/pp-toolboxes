function [cfg,skipprocessing] = lab_set_reference_data(cfg,header)

disp ('   Ask for new-references-settings')

skipprocessing = 0;

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('cfg','var') | ~isfield(cfg,'REF') | isempty(cfg.REF)
    cfg.REF.folder = 'REF';
    cfg.REF.format = 'sef';
    cfg.REF.eegsource = {'mean'};
    cfg.REF.interpolatebad=true;
    cfg.REF.LAPL = [];
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

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Reference','eegsource'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = refitems;
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [80 90];
Formats(end,1).callback = {@lab_get_eegsource,'@ALL','@ALL',cfg,header};
Formats(end,1).span = [2 1];

Prompt(end+1,:) = {'Montage','montage'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_load_montage,'montage','montage',cfg,header};

Prompt(end+1,:) = {'Laplacian','LAPL'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_laplacian,'LAPL','LAPL'};

Prompt(end+1,:) = {'Interpolate bad' 'interpolatebad'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'File format','format'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input';
Formats(end,1).items = {'sef';'edf';'eph';'ep';'txt'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [60 85];
Formats(end,1).callback = {@lab_get_format,'scaletxt','format','scaletxt'};

Prompt(end+1,:) = {'Scale TXT','scaletxt'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_show_result,'scaletxt','scaletxt'};

[cfg.REF,Cancelled] = inputsdlg(Prompt,'Reference data',Formats,cfg.REF);
if isempty(cfg.REF) | Cancelled == 1
    cfg.REF = [];
    skipprocessing = 1;
    return
end