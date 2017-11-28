function settings = lab_set_plot_ISresult

settings.mri_file = '';
settings.ris_file = '';
settings.sizedots = 4;
settings.dotcolor = [1 0 0];

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'MRI-file','mri_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr;*.nii','MRI-file (.hdr/.nii)';'*.*','All Files'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = 250;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'IS-Result-File','ris_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.ris','Select IS-Result-File'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = 250;

Prompt(end+1,:) = {'Size','sizedots'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).size = 40;
Formats(end,1).limits = [0 inf];

Prompt(end+1,:) = {'Color','dotcolor'};
Formats(end+1,1).type = 'color';

[settings,Cancelled] = inputsdlg(Prompt,'Plot IS results',Formats,settings);
if Cancelled == 1
    settings = [];
    return
end

end