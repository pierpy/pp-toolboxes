function [cfg,skipprocessing] = lab_set_detect_dbs(cfg,calc)

skipprocessing = 0;

if ~exist('cfg','var')
    cfg = [];
end

if exist('calc','var') & isfield(calc,'Filelist')
    Output = lab_prepare_subjectname(calc.Filelist{1,1});
else
    Output = {[],[]};
end
if ~isfield(cfg,'DBS') | ~isfield(cfg.DBS,'subjectname') | isempty(cfg.DBS.subjectname)
    if ~isempty(Output{1})
        cfg.DBS.subjectname = -1;
    else
        cfg.DBS.subjectname = 0;
    end
end
if ~isfield(cfg.DBS,'folder') | isempty(cfg.DBS.folder)
    cfg.DBS.folder = 'DBS';
    cfg.DBS.threshold = 2000;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Number of underscores in subject name',''};
Formats(end+1,1).type = 'text';

if ~isempty(Output{1})
    Prompt(end+1,:) = {[Output{1} ' ' Output{2}],'subjectname'};
else
    Prompt(end+1,:) = {Output{2},'subjectname'};
end
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-99 99];
Formats(end,1).size = 40;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 100;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Power - Threshold','threshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 50;

[cfg.DBS,Cancelled] = inputsdlg(Prompt,'Detect DBS',Formats,cfg.DBS);
pause(0.1);
if isempty(cfg.DBS) | Cancelled == 1
    cfg.DBS = [];
    skipprocessing = 1;
end

end