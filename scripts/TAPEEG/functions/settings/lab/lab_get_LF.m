function LF = lab_get_LF(LF,IS_file,LOC_file)

if ~exist('LOC_file','var')
    LOC_file = '';
end
if ~exist('IS_file','var')
    IS_file = '';
end

dodialog = false;
switch IS_file
    case 'MRI-file (iso.fif)'
        ISformat = 'fif';
        dodialog = true;
    case 'MRI-file (mri.fif)'
        ISformat = 'fif';
        dodialog = true;
    case 'MRI-file (.hdr)'
        ISformat = 'hdr';
        dodialog = true;
    case 'MRI-file (.nii)'
        ISformat = 'nii';
        dodialog = true;
    case 'MRI-file (dicom)'
        ISformat = 'dicom';
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
    LF = [];
    return
end

if ~exist('LF','var') | ~isfield(LF,'method')
    if strcmp(ISformat,'fif') | strcmp(ISformat,'fiff')
        LF.method = 'localspheres';
    else
        LF.method = 'openmeeg';
    end
end
if ~isfield(LF,'skullconduct')
    if ~strcmp(LF.method,'FEM-simbio')
        LF.skullconduct = 0.033;
    else
        LF.skullconduct = [0.33 0.14 1.79 0.01 0.431];
    end
end
if ~isfield(LF,'coregmode')
    LF.coregmode = 0;
end
if ~isfield(LF,'correctdigits')
    LF.correctdigits = false;
end
if ~isfield(LF,'interactive')
    LF.interactive = false;
end
if ~isfield(LF,'maxdist')
    LF.maxdist = 1;
end

if ~isfield(LF,'MESH') | isempty(LF.MESH)
    LF.MESH.Nscalp = 3000;
    LF.MESH.Nskull = 3000;
    LF.MESH.Nbrain = 3000;
    LF.MESH.SkullDist = 4;
    LF.MESH.SkullNr = 0;
    LF.MESH.OtherDist = 1;
    LF.MESH.shift = 0.3;
    LF.MESH.resolution = 1;
end

if isfield(LF,'coregmode') & ~isempty(LF.coregmode)
    LF.coregmode = LF.coregmode + 1;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Method for Leadfield-calculation','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'singlesphere','localspheres','openmeeg','dipoli','concentricspheres','FEM-simbio'};
Formats(end,1).callback = {@control_method,'@ALL','@ALL'};

Prompt(end+1,:) = {'settings','MESH'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_MESH,'MESH','MESH','method'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Skull conductivity','skullconduct'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 250;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Coregister electrodes (EEG)' 'coregmode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'Off','Automatic rotation/translation','Interactive coregistration','Use landmarks'};
Formats(end,1).callback = {@set_landmarks,'landmarks','landmarks',LOC_file,'coregmode',IS_file};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Landmarks' 'landmarks'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_landmarks,'landmarks','landmarks',LOC_file,'coregmode',IS_file};

Prompt(end+1,:) = {'Eliminate electrodes after coregistration' 'deletelocs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@select_electrodes,'deletelocs','deletelocs',LOC_file};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Correct Digitizer-points (MEG)' 'correctdigits'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

[LF,Cancelled] = inputsdlg(Prompt,'LF settings',Formats,LF);
if Cancelled == 1
    LF = [];
else
    LF.coregmode = LF.coregmode - 1;
end

end

function settings = control_method(settings)
   if strcmp(settings.method,'dipoli')
       settings.MESH.Nscalp = 3000;
       settings.MESH.Nskull = 3000;
       settings.MESH.Nbrain = 5000;
   else
       settings.MESH.Nscalp = 3000;
       settings.MESH.Nskull = 3000;
       settings.MESH.Nbrain = 3000;
   end
   if strcmp(settings.method,'FEM-simbio') & length(settings.skullconduct) ~= 5
       settings.skullconduct = [0.33 0.14 1.79 0.01 0.431];
   elseif ~strcmp(settings.method,'FEM-simbio') & length(settings.skullconduct) ~= 1
       settings.skullconduct = 0.033;
   end
   settings.MESH = lab_get_MESH(settings.MESH,settings.method);
end

function selection = select_electrodes(selection,LOC_file)
    if ischar(LOC_file) & exist(LOC_file,'file')
        settings.LOCS = lab_read_locs(LOC_file);
    else
        settings.LOCS = lab_read_locs;
    end
    if isempty(settings.LOCS)
        selection = [];
        return
    end
    settings.indexed = selection;
    settings.Color = [1 1 1];
    settings.ColorIdx = [1 0 0];
    settings.Title = 'Select Channels';
    selection = lab_plot_locs(settings,1);
end

function settings = set_landmarks(settings,LOC_file,coregmode,MRI_file)
    if exist('coregmode','var') & coregmode ~= 4
        settings = [];
        return
    end
    if ~isfield(settings,'landmarks')
        settings.landmarks = [];
        settings.mrilandmarks = [];
    end
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Select Electrodes','landmarks'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).callback = {@select_electrodes,'landmarks','landmarks',LOC_file};
    Formats(end,1).enable = 'inactive';
    
    Formats(end+1,1).type = 'none';
    
    if exist('MRI_file','var') & exist(MRI_file,'file')
        if isempty(settings.mrilandmarks)
            settings.mrilandmarks = MRI_file;
        end
        Prompt(end+1,:) = {'Edit landmarks', 'mrilandmarks'};
    else
        Prompt(end+1,:) = {'Template-MRI with landmarks', 'mrilandmarks'};
    end
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'text';
    Formats(end,1).size = 200;
    Formats(end,1).enable = 'inactive';
    Formats(end,1).callback = {@edit_landmarks,'mrilandmarks','mrilandmarks',LOC_file,'landmarks'};
    
    [settings,Cancelled] = inputsdlg(Prompt,'Set landmarks',Formats,settings);
    if Cancelled == 1
        settings = [];
    end
    return
end

function mrifile = edit_landmarks(mrifile,LOC_file,landmarks)
    if isempty(mrifile) | ~exist(mrifile,'file')
        [tmp1,tmp2] = uigetfile('*.hdr','Select MRI-Template-file with landmarks');
        mrifile = fullfile(tmp2,tmp1);
        clearvars tmp1 tmp2
    elseif isempty(landmarks)
        button = questdlg('Delete Template-MRI','Delete Template-MRI for definition of landmarks','Yes','No','No');
        if strcmp(button,'Yes')
            mrifile = '';
        end
        return
    end
    settings.landmarks = landmarks;
    settings.forceselection = true;
    if ischar(LOC_file) & exist(LOC_file,'file')
        LOCS = lab_read_locs(LOC_file);
    else
        LOCS = lab_read_locs;
    end
    mrifile = lab_mri_landmarks(mrifile,LOCS,settings);
end