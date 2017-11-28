function settings = lab_set_plot_mri(settings,flag)

if ~exist('flag','var')
    flag = false;
end
if ~exist('settings','var') | ~isfield(settings,'MRI_file')
    settings.MRI_file = '';
    settings.ATLAS_file = '';
    settings.SPI_file = '';
    settings.RIS_file = '';
    settings.Interpolate = 0;
    settings.color = [1 0 0];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'MRI-file','MRI_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii;*.fif;*.fiff;*.img;*.dcm','MRI-file'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];

if flag == true
    Formats(end+1,1).type = 'none';
end

Prompt(end+1,:) = {'Atlas-File','ATLAS_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii','MRI-file (.hdr/.nii)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];

if flag == false
    Formats(end+1,1).type = 'none';
    
    Prompt(end+1,:) = {'Solutionpoints-file','SPI_file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.spi','Solutionpoints (.spi)'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).size = [400 0];
    
    Prompt(end+1,:) = {'Result-file','RIS_file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.ris','Inverse solution results (.ris)'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).size = [400 0];
    
    Prompt(end+1,:) = {'Interplate factor (0 = auto)','Interpolate'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {0,1,2,3,4,5,6,7,8,9,10,15,20};
    
    Prompt(end+1,:) = {'Color','color'};
    Formats(end+1,1).type = 'color';
end

[settings,Cancelled] = inputsdlg(Prompt,'Plot MRI',Formats,settings);
if Cancelled == 1
    settings = [];
end