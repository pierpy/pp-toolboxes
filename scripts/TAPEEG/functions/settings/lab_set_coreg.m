function [settings,skipprocessing] = lab_set_coreg(settings,header,skiploc,domrifile)

skipprocessing = 0;

if ~exist('domrifile','var')
    domrifile = true;
end
if ~exist('skiploc','var')
    skiploc = false;
end

global THeader
if ~isempty(THeader)
    header = THeader;
elseif ~exist('header','var')
    header = [];
end

if ~exist('settings','var') | ~isfield(settings,'coregmode')
    settings.coregmode = 1;
    settings.deletelocs = [];
    settings.landmarks = [];
    settings.maxdist = 1;
    settings.correctdigits = false;
    settings.findsphere = false;
end

if isfield(settings,'MRI_file')
    switch settings.MRI_file
        case '.hdr'
            settings.MRI_file = 'MRI-file (.hdr)';
        case '.nii'
            settings.MRI_file = 'MRI-file (.nii)';
        case 'iso.fif'
            settings.MRI_file = 'MRI-file (iso.fif)';
        case 'mri.fif'
            settings.MRI_file = 'MRI-file (mri.fif)';
        case 'dicom'
            settings.MRI_file = 'MRI-file (dicom)';
    end
else
    settings.MRI_file = '';
end

Prompt = cell(0,2);
Formats = {};

Prompt{end+1,1} = 'MRI-File';
Formats(end+1,1).type = 'text';
Formats(end,1).style = 'text';
Formats(end,1).span = [1 3];

if domrifile == false
    Prompt(end+1,:) = {'','MRI_file'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'','Select File','MRI-file (.hdr)','MRI-file (.nii)', ...
        'MRI-file (iso.fif)','MRI-file (mri.fif)','MRI-file (dicom)'};
    Formats(end,1).size = 300;
    Formats(end,1).callback = {@get_mrifile,'@ALL','@ALL'};
    Formats(end,1).span = [1 3];
else
    Prompt(end+1,:) = {'','MRI_file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.hdr;*.nii;*.fif;*.fiff','Select MRI-file'};
    Formats(end,1).limits = [0 1];
    Formats(end,1).size = [300 0];
    Formats(end,1).span = [1 3];
end

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

if skiploc == false
    Prompt{end+1,1} = 'Locations-File';
    Formats(end+1,1).type = 'text';
    Formats(end,1).style = 'text';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'','LOC_file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.xyz;*.els','Electrodes-file (.els/.xyz)'; ...
        '*.fif;*.fiff','Data-file (.fif)'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).size = [300 0];
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
else    
    settings.LOC_file = '';
end

Prompt(end+1,:) = {'Mode','coregmode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'Automatic rotation/translation','Interactive coregistration','Use landmarks'};
Formats(end,1).callback = {@select_electrodes,'landmarks','landmarks','LOC_file',header,'coregmode','MRI_file'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Landmarks' 'landmarks'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@select_electrodes,'landmarks','landmarks','LOC_file',header,'coregmode','MRI_file'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Eliminate electrodes after coregistration' 'deletelocs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@select_electrodes,'deletelocs','deletelocs','LOC_file',header};
Formats(end,1).span = [1 3];

if isempty(header) | isfield(header.locs,'digits')
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'Correct Digitizer-points' 'correctdigits'};
    Formats(end+1,1).type = 'check';
        
    Prompt(end+1,:) = {'Maximal distance for nose points','maxdist'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [0 inf];
    Formats(end,1).size = 30;
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Find optimal sphere' 'findsphere'};
    Formats(end+1,1).type = 'check';
end

[settings,Cancelled] = inputsdlg(Prompt,'Coregistration',Formats,settings);
if Cancelled == 1 | isempty(settings)
    skipprocessing = 1;
    settings = [];
else
    switch settings.MRI_file
        case 'MRI-file (.hdr)'
            settings.MRI_file = '.hdr';
        case 'MRI-file (.nii)'
            settings.MRI_file = '.nii';
        case 'MRI-file (iso.fif)'
            settings.MRI_file = 'iso.fif';
        case 'MRI-file (mri.fif)'
            settings.MRI_file = 'mri.fif';
        case 'MRI-file (dicom)'
            settings.MRI_file = 'dicom';
    end
    if ~isempty(settings.LOC_file) & length(settings.LOC_file) > 2 & strcmp(settings.LOC_file(end-2:end),'fif')
        settings.correctdigits = true;
        settings.findsphere = true;
    end
end

end

function settings = get_mrifile(settings)
    if strcmp(settings.MRI_file,'Select File')
        [MRI_file,MRI_filepath] = uigetfile({'*.hdr;*.nii;*.fif;*.fiff','MRI-file (.hdr/.nii/.fif)'; ...
            '*.img;*.ima;*.dcm','DICOM (.img/.ima/.dcm)'},'Select MRI_file');
        if ischar(MRI_file) & exist(fullfile(MRI_filepath,MRI_file),'file')
            settings.MRI_file = fullfile(MRI_filepath,MRI_file);
        else
            settings.MRI_file = '';
        end
    end
end

function selection = select_electrodes(selection,LOC_file,header,coregmode,MRI_file)
    if exist('coregmode','var') & coregmode ~= 3
        selection = [];
        return
    end
    if ischar(LOC_file) & exist(LOC_file,'file')
        settings.LOCS = lab_read_locs(LOC_file);
    elseif isfield(header,'locs') & ~isempty(header.locs)
        settings.LOC = header.locs;
    else
        settings.selection = selection;
        Prompt = {'Electrodes','selection'};
        Formats.type = 'edit';
        Formats.format = 'vector';
        Formats.limits = [-inf inf];
        [settings,Cancelled] = inputsdlg(Prompt,'Select electrodes',Formats,settings);
        if Cancelled == 0
            selection = settings.selection;
        else
            selection = [];
        end
        return
    end
    settings.indexed = selection;
    settings.Color = [1 1 1];
    settings.ColorIdx = [1 0 0];
    settings.Title = 'Select Channels';
    selection = lab_plot_locs(settings,1);
    if exist('MRI_file','var') & exist(MRI_file,'file')
        settings.landmarks = selection;
        settings.forceselection = true;
        lab_mri_landmarks(MRI_file,settings.LOCS,settings);
    end
end