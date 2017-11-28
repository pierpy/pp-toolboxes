function [settings,skipprocessing] = lab_set_match_atlas2mri(settings)

disp ('   Ask for settings')
skipprocessing = 0;

if ~exist('settings','var')
    settings = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Atlas-file','Atlas_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii','Atlas-file (.hdr/.nii)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'MRI-file','MRI_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii','MRI-file (.hdr/.nii)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Atlas-Template','AtlasTemplate_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii','AtlasTemplate-file (.hdr/.nii)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];

[settings,Cancelled] = inputsdlg(Prompt,'Match Atlas2MRI',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
end

