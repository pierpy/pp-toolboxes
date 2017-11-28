function [settings,skipprocessing] = lab_set_plot_mesh(flagbnd)

if ~exist('flagbnd','var')
    flagbnd = true;
end

settings.plotfaces = true;
settings.facecolor = [0 1 0];
settings.alpha = 0.5;
settings.gapfactor = 2;
settings.plotedges = false;
settings.edgecolor = [0 0 0];
settings.plotdots = false;
settings.dotscolor = [0 0 1];
settings.sizedots = 4;
settings.plotlabels = false;
settings.labelsize = 20;
settings.labeldistance = 1;

Prompt = {};
Formats = {};

if flagbnd == true
    settings.meshfile = [];
    Prompt(end+1,:) = {'Mesh file','meshfile'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.els;*.xyz;*.grad','Electrodes-file (.els/.xyz)'; ...
        '*.spi','Solutionpoints (.spi)';'*.hdr;*.nii','MRI-file (.hdr/.nii)';'*.*','All Files'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).size = [300 0];
    Formats(end,1).span = [1 4];

    Prompt(end+1,:) = {'',''};
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 4];
end

Prompt(end+1,:) = {'Plot faces','plotfaces'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Color','facecolor'};
Formats(end+1,1).type = 'color';
Formats(end,1).format = 'color';

Prompt(end+1,:) = {'Alpha','alpha'};
Formats(end+1,1).type = 'range';
Formats(end,1).limits = [0 1];

Prompt(end+1,:) = {'Gap factor', 'gapfactor'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Plot edges','plotedges'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Color','ecolor'};
Formats(end+1,1).type = 'color';
Formats(end,1).format = 'color';

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Plot dots','plotdots'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Color','dotscolor'};
Formats(end+1,1).type = 'color';
Formats(end,1).format = 'color';

Prompt(end+1,:) = {'Size', 'sizedots'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Plot labels','plotlabels'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Size', 'labelsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Distance', 'labeldistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

[settings,skipprocessing] = inputsdlg(Prompt,'Plot Mesh',Formats,settings);
pause(0.2);

if skipprocessing ~= 1
    if isfield(settings,'meshfile') & exist(settings.meshfile,'file')
        settings.bnd = lab_read_data(settings.meshfile);
        if isfield(settings.bnd,'img') | isfield(settings.bnd,'anatomy')
            settings.bnd = lab_mesh_vol(settings.bnd,10000);
        end
    end
end