function settings = lab_get_ROIS(settings)

if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'SPI_file')
    settings.SPI_file = '';
    SPIformat = '';
elseif ~isempty(settings.SPI_file)
    switch settings.SPI_file
        case 'Solutionpoints (.spi)'
            SPIformat = 'spi';
        case 'MRI-file (IS-file)'
            SPIformat = 'hdr';
        case 'MRI-file (hdr)'
            SPIformat = 'hdr';
        case 'MRI-file (nii)'
            SPIformat = 'nii';
        otherwise
            [~,~,SPIformat] = lab_filename(settings.SPI_file);
    end
else
    SPIformat = '';
end

if ~isfield(settings,'ROIS_file') | isempty(settings.ROIS_file) | strcmp(settings.ROIS_file,'Select File')
    if any(strcmp(settings.ROISmethods,'pseudoelectrodes'))
        [ROIS_file,ROIS_filepath] = uigetfile({'*.xyz;*.els','LOC-file (.xyz/.els)'},'Select Electrodes-file');
    elseif strcmp(SPIformat,'hdr') | strcmp(SPIformat,'nii')
        [ROIS_file,ROIS_filepath] = uigetfile({'*.hdr;*.nii','MRI-file (.hdr/.nii)'},'Select Atlas-file');
    elseif strcmp(SPIformat,'spi')
        [ROIS_file,ROIS_filepath] = uigetfile({'*.rois','ROIS-file (.rois)'}, ...%edit G Bogaarts 2016
            'Select ROIS-file');
    else
        [ROIS_file,ROIS_filepath] = uigetfile({'*.rois','ROIS-file (.rois)'; ...
            '*.hdr;*.nii','MRI-file (.hdr/.nii)';'*.xyz;*.els','LOC-file (.xyz/.els)'},'Select ROIS-file');
    end
    settings.ROIS_file = fullfile(ROIS_filepath,ROIS_file);
end

switch settings.ROIS_file
    case 'MRI-Atlas (.nii)'
        ROISformat = 'nii';
    case 'MRI-Atlas (.hdr)'
        ROISformat = 'hdr';
    case 'ROIS-file (.rois)'
        ROISformat = 'rois';
    otherwise
        [~,~,ROISformat] = lab_filename(settings.ROIS_file);
end

if (strcmp(ROISformat,'hdr') | strcmp(ROISformat,'nii')) & strcmp(SPIformat,'spi')
    if strcmp(settings.SPI_file,'Solutionpoints (.spi)')
        settings.ROIS_file = 'ROIS-file (.rois)';
    else
        settings.ROIS_file = '';
    end
    ROISformat = 'rois';
elseif strcmp(ROISformat,'rois') & (strcmp(SPIformat,'hdr') | strcmp(SPIformat,'nii'))
    settings.ROIS_file = '';
    settings.ROISsettings = [];
    return
end

if ~isfield(settings,'ROISsettings')
    settings.ROISsettings = [];
end

if strcmp(ROISformat,'hdr') | strcmp(ROISformat,'nii')
    Prompt = {};
    Formats = {};
    if isempty(settings.ROISsettings)
        settings.ROISsettings.AtlasTemplate = '';
        settings.ROISsettings.CoregAtlas = false;
        settings.ROISsettings.ReduceAAL = true;
        settings.ROISsettings.RemoveBorder = true;
    end
    if strcmp(SPIformat,'hdr') | strcmp(SPIformat,'nii') | strcmp(SPIformat,'fif') | strcmp(SPIformat,'fiff')
        settings.ROISsettings.CoregAtlas = true;
        if exist(settings.ROIS_file,'file')
            [~,ROIS_filepath] = lab_filename(settings.ROIS_file);
            if exist(fullfile(ROIS_filepath,'Brain.nii'),'file')
                settings.ROISsettings.AtlasTemplate = fullfile(ROIS_filepath,'Brain.nii');
            end
        end
        
        Prompt(end+1,:) = {'Atlas Brain-Template','AtlasTemplate'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).callback = {@get_template,'AtlasTemplate','AtlasTemplate'};
        
        Prompt(end+1,:) = {'Coregister Atlas to MRI-file','CoregAtlas'};
        Formats(end+1,1).type = 'check';
        settings.CoregAtlas = true;
        
        Formats(end+1,1).type = 'none';
        Formats(end,1).span = [1 2];
    end
    Prompt(end+1,:) = {'Remove borders of regions','RemoveBorder'};
    Formats(end+1,1).type = 'check';
    
    Prompt(end+1,:) = {'Reduce AAL-Atlas to standard regions','ReduceAAL'};
    Formats(end+1,1).type = 'check';
    
    settings.ROISsettings = inputsdlg(Prompt,'ROIS Extra',Formats,settings.ROISsettings,2);
else
    settings.ROISsettings = [];
end

end

function Template = get_template(Template)
    settings.Template = Template;
    Prompt = {'Template for normalization','Template'};
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