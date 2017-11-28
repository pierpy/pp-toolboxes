function [settings,skipprocessing] = lab_set_compute_leadfield(settings)
disp ('   Ask for compute LF settings')
skipprocessing = 0;

if ~exist('settings','var')
    settings = [];
    if ~isfield(settings,'method')
        settings.SPI.method = 'Gray';
    end
    if ~isfield(settings,'numspi')
        settings.SPI.numspi = 10000;
    end
    if ~isfield(settings,'atlasfile')
        settings.SPI.atlasfile = '';
    end
    if ~isfield(settings,'RemoveBorder')
        settings.SPI.RemoveBorder = false;
    end
    if ~isfield(settings,'atlasexclude')
        settings.SPI.atlasexclude = [];
    end
    if ~isfield(settings,'ReduceAAL')
        settings.SPI.ReduceAAL = false;
    end
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'MRI-file','mrifile'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr','MRI-file (.hdr)';'*.mat','Headmodel (.mat)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [400 0];
Formats(end,1).callback = {@lab_get_LF,'LF','LF','mrifile','LOC_file'};

Prompt(end+1,:) = {'MRI settings','MRI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_pMRI,{'MRI','LF'},'MRI','mrifile','LF'};

Prompt(end+1,:) = {'Solutionpoints-file','SPI_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr','MRI-file (.hdr)';'*.spi','Solutionpoints (.spi)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [-1 0];
Formats(end,1).callback = {@lab_get_spi,'SPI','SPI','SPI_file','mrifile'};

Prompt(end+1,:) = {'SPI-settings','SPI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_spi,'SPI','SPI','SPI_file','mrifile'};

Prompt(end+1,:) = {'Electrodes-file','LOC_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.els;*.xyz','Electrodes-file (.els/.xyz)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [-1 0];

Prompt(end+1,:) = {'Leadfield','LF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_LF,'LF','LF','mrifile','LOC_file'};
Formats(end,1).span = [1 2];

[settings,Cancelled] = inputsdlg(Prompt,'Compute Leadfield',Formats,settings,2);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
elseif isfield(settings,'LF') & ~isempty(settings.LF) & isfield(settings.LF,'skullconduct')
    settings.skullconduct = settings.LF.skullconduct;
    settings.coregmode = settings.LF.coregmode;
    settings.landmarks = settings.LF.landmarks;
    settings.correctdigits = settings.LF.correctdigits;
    settings.deletelocs = settings.LF.deletelocs;
    settings.maxdist = settings.LF.maxdist;
    settings.method = settings.LF.method;
    if isfield(settings.LF,'MESH')
        settings.MESH = settings.LF.MESH;
    end
else
    settings.LF = [];
end

return