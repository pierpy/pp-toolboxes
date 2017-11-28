function [cfg,skipprocessing] = lab_set_process_mri(cfg,setsearch)

skipprocessing = 0;
if ~exist('setsearch','var')
    setsearch = false;
end

if ~exist('cfg','var') | ~isfield(cfg,'MRI')
    cfg.MRI.FILT = [];
    cfg.MRI.template = '';
    cfg.MRI.matchmode = 'Coregister';
    cfg.MRI.brainthreshold = 0.5;
    cfg.MRI.scalpthreshold = 0.1;
    cfg.MRI.SEGbrain = true;
    cfg.MRI.SEGcorrect = false;
    cfg.MRI.atlas = '';
    cfg.MRI.collectsegm = 'gray';
    cfg.MRI.collectthreshold = 0.5;
    cfg.MRI.SPI = [];
    cfg.MRI.LOC_file = '';
    cfg.MRI.LF = [];
    cfg.MRI.SEARCH.searchstring{1,1} = '.hdr';
end
if ~isfield(cfg.MRI,'orientation') | isempty(cfg.MRI.orientation) | length(cfg.MRI.orientation) < 3
    cfg.MRI.orientation = [1 2 3];
    cfg.MRI.Orient = 'RAS';
else
    Oletters = 'RASLPI';
    cfg.MRI.Orient = Oletters(cfg.MRI.orientation);
end
if isfield(cfg,'EEG_file') & isfield(cfg,'EEG_filepath')
    MRI_file = fullfile(cfg.EEG_filepath,cfg.EEG_file);
else
    MRI_file = '';
end

Prompt = cell(0,2);
Formats = [];

if setsearch == true
    Prompt(end+1,:) = {'Search Files','SEARCH'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_set_searchstrings,'@ALL','@ALL',0,0,{'strings','MRI (.hdr)'}};
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
end

Prompt(end+1,:) = {'Input orientation','Orient'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [30 20];
Formats(end,1).callback = {@set_orientation,'Orient','Orient',MRI_file};

Prompt(end+1,:) = {'Set contrast','CONTRAST'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_contrast,'CONTRAST','CONTRAST',MRI_file};

Prompt(end+1,:) = {'Anisotropic filtering','FILT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_set_mrifiltering,'FILT','FILT'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Match to Template','template'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii;*.fif','Select template'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = 200;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Mode','matchmode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Coregister','Normalize'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Segment MRI',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'gray/white/csf','SEGbrain'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Threshold brain', 'brainthreshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 1];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Threshold scalp', 'scalpthreshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 1];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Correct segmentation','SEGcorrect'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Probability maps','SEGprobmaps'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Calculate leadfield',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'SPI-settings','SPI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_spi,'SPI','SPI','MRI-file (.hdr)','MRI-file (.hdr)'};

Prompt(end+1,:) = {'Electrodes-file','LOC_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.els;*.xyz','Electrodes-file (.els/.xyz)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = 200;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Leadfield','LF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_LF,'LF','LF','MRI-file (.hdr)','LOC_file'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Collect Voxels',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'', 'collectsegm'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'','gray','white','csf','brain','skull','scalp'};

Prompt(end+1,:) = {'Threshold', 'collectthreshold'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 1];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Atlas','atlas'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii;*.fif','Select atlas'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = 250;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Atlas Brain-Template','AtlasTemplate'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_template,'AtlasTemplate','AtlasTemplate','atlas'};

[cfg.MRI,Cancelled] = inputsdlg(Prompt,'Process mri',Formats,cfg.MRI);
if Cancelled == 1
    skipprocessing = 1;
    return
else
    Oletters = 'RASLPI';
    for i = 1:3
        tmp = strfind(Oletters,cfg.MRI.Orient(i));
        if ~isempty(tmp)
            cfg.MRI.orientation(1,i) = tmp(1);
        else
            cfg.MRI.orientation(1,i) = 0;
        end
    end
    if setsearch == true & ~isempty(cfg.MRI.SEARCH)
        cfg.SEARCH = cfg.MRI.SEARCH;
    end
end

end

function Orient = set_orientation(Orient,MRI_file)

Oletters = 'RASLPI';
if length(Orient) == 3
    for i = 1:3
        tmp = strfind(Oletters,Orient(i));
        if ~isempty(tmp)
            orientation(1,i) = tmp(1);
        else
            orientation(1,i) = 0;
        end
    end
else
    orientation = [1 2 3];
end

if exist(MRI_file,'file')
    mri = lab_read_mri(MRI_file);
else
    mri = lab_read_mri;
end

if ~isfield(mri,'anatomy')
    orientation = lab_get_orient(orientation);
else
    [mri,settings] = lab_correct_mri2RAS(mri);
    if ~isempty(mri)
        orientation = settings.orient;
    end
end
Orient = Oletters(orientation);

end

function contrast = set_contrast(contrast,MRI_file)
    if ~isempty(MRI_file) & exist(MRI_file,'file')
        mri = lab_read_mri(MRI_file);
    else
        mri = [];
    end
    if isfield(mri,'anatomy')
        settings2.fhandle = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off', ...
            'Menubar','none','Name','Correct Orientation to RAS');
        pos = get(0,'ScreenSize');
        pos = [100 (pos(4)-800) 700 700];
        if pos(2) < 0
            pos(2) = 100;
        end
        set(settings2.fhandle,'Position',pos);
        settings2.setminmax = true;
        if isfield(contrast,'minval')
            mri.minval = contrast.minval;
        end
        if isfield(contrast,'maxval')
            mri.maxval = contrast.maxval;
        end
        mri = lab_plot_orthoslides(mri,settings2);
        close(settings2.fhandle);
        if isfield(mri,'minval') & isfield(mri,'maxval')
            contrast.minval = mri.minval;
            contrast.maxval = mri.maxval;
        end
    else
        Prompt = {'Minimal Value','minval';'Maximal Value','maxval'};
        Formats.type = 'edit';
        Formats.format = 'float';
        Formats.limits = [0 inf];
        Formats.size = 40;
        Formats.span = [1 2];
        Formats = [Formats;Formats];
        [contrast,Cancelled] = inputsdlg(Prompt,'Set Min/Max Values',Formats,contrast);
        if Cancelled == 1
            contrast = [];
        end
    end
end

function Template = get_template(Template,atlas)
    settings.Template = Template;
    if ~isempty(atlas) & isempty(Template)
        [~,atlas_filepath] = lab_filename(atlas);
        if exist(fullfile(atlas_filepath,'Brain.nii'),'file')
            settings.Template = fullfile(atlas_filepath,'Brain.nii');
        end
    end
    Prompt = {'Template for atlas normalization','Template'};
    Formats.type = 'edit';
    Formats.format = 'file';
    Formats.items = {'*.nii','Brain-MRI-file (*.nii)'};
    Formats.limits = [0 1];
    Formats.size = 300;
    [settings,Cancelled] = inputsdlg(Prompt,'Normalization Template',Formats,settings);
    if Cancelled == 1
        Template = '';
    else
        Template = settings.Template;
    end
end