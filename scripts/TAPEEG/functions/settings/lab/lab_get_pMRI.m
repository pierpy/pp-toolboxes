function [MRI,LF] = lab_get_pMRI(MRI,IS_file,LF,IS_files)

if ~exist('LF','var')
    LF = [];
end
if ~exist('IS_file','var')
    IS_file = '';
end

dodialog = false;
switch IS_file
    case 'MRI-file (iso.fif)'
        dodialog = true;
    case 'MRI-file (mri.fif)'
        dodialog = true;
    case 'MRI-file (.hdr)'
        dodialog = true;
    case 'MRI-file (.nii)'
        dodialog = true;
    case 'MRI-file (dicom)'
        dodialog = true;
    otherwise
        [~,~,ISformat] = lab_filename(IS_file);
        if isempty(ISformat)
            ISformat = IS_file;
        end
        if strcmp(ISformat,'fif') | strcmp(ISformat,'fiff') | strcmp(ISformat,'hdr') | strcmp(ISformat,'nii') | strcmp(ISformat,'dicom')
            dodialog = true;
        end
end

if dodialog == false
    MRI = [];
    return
end

if exist(IS_file,'file')
    MRI_file = IS_file;
elseif exist('IS_files','var') & isfield(IS_files,'MRI_file') & exist(IS_files.MRI_file,'file')
    MRI_file = IS_files.MRI_file;
else
    MRI_file = '';
end
    

if ~exist('MRI','var') | ~isfield(MRI,'orientation') | isempty(MRI.orientation) | length(MRI.orientation) < 3
    MRI.orientation = [1 2 3];
    MRI.Orient = 'RAS';
    MRI.FILT = [];
else
    tmp = 'RASLPI';
    MRI.Orient = tmp(MRI.orientation);
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Input orientation','Orient'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [80 20];
Formats(end,1).callback = {@set_orientation,'Orient','Orient',MRI_file};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Set contrast','CONTRAST'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_contrast,'CONTRAST','CONTRAST',MRI_file};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Anisotropic filtering','FILT'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_mrifiltering,'FILT','FILT'};

[MRI,Cancelled] = inputsdlg(Prompt,'MRI settings',Formats,MRI);
if Cancelled == 1
    MRI = [];
    if isfield(LF,'MRI')
        LF.MRI = [];
    end
else
    tmp = 'RASLPI';
    for i = 1:3
        tmp2 = strfind(tmp,MRI.Orient(i));
        if ~isempty(tmp2)
            MRI.orientation(1,i) = tmp2(1);
        else
            MRI.orientation(1,i) = 0;
        end
    end
    LF.MRI = MRI;
end

end

function settings = set_mrifiltering(settings)

if isempty(settings)
    settings.ts = 0.0625;
    settings.iter = 5;
    settings.cond = 3;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Sampling rate', 'ts'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 10];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Number of iterations', 'iter'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Diffusion parameter', 'cond'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 30;

[settings,Cancelled] = inputsdlg(Prompt,'MRI filter',Formats,settings);
if Cancelled == 1
    settings = [];
    return
end

end

function Orient = set_orientation(Orient,MRI_file)

tmp = 'RASLPI';
if length(Orient) == 3
    for i = 1:3
        tmp2 = strfind(tmp,Orient(i));
        if ~isempty(tmp2)
            orientation(1,i) = tmp2(1);
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