function [settings,SPI_file] = lab_get_spi(settings,SPI_file,IS_file)

if exist('IS_file','var')
    switch IS_file
        case 'Cartool (.is)'
            ISformat = 'is';
        case 'sLoreta (.spinv)'
            ISformat = 'spinv';
        case 'Leadfield (LF.bin)'
            ISformat = 'bin';
        case 'Headmodel (.mat)'
            ISformat = 'mat';
        case 'MRI-file (.hdr)'
            ISformat = 'hdr';
        case 'MRI-file (.nii)'
            ISformat = 'nii';
        case 'MRI-file (iso.fif)'
            ISformat = 'fif';
        case 'MRI-file (mri.fif)'
            ISformat = 'fif';
        case 'MRI-file (dicom)'
            ISformat = 'dicom';
        otherwise
            [~,~,ISformat] = lab_filename(IS_file);
    end
else
    ISformat = '';
end

if exist('SPI_file','var') & ~isempty(SPI_file)
    if strcmp(SPI_file,'Select File')
        if ~any(strcmp({'is','spinv','bin','mat'},ISformat))
            [SPI_file,SPI_filepath] = uigetfile({'*.hdr;*.nii','MRI-file (.hdr/.nii/.fif)'; ...
                '*.img;*.ima;*.dcm','DICOM (.img/.ima/.dcm)'; ...
                '*.spi','Solutionpoints (.spi)'},'Select SPI-file');
        else
            [SPI_file,SPI_filepath] = uigetfile({'*.spi','Solutionpoints (.spi)'}, ...
                'Select SPI-file');
        end
        SPI_file = fullfile(SPI_filepath,SPI_file);
    end
else
    SPI_file = '';
    settings = [];
    return
end

switch SPI_file
    case 'Solutionpoints (.spi)'
        SPIformat = 'spi';
    case 'MRI-file (IS-file)'
        SPIformat = 'hdr';
    case 'MRI-file (hdr)'
        SPIformat = 'hdr';
    case 'MRI-file (nii)'
        SPIformat = 'nii';
    otherwise
        [~,~,SPIformat] = lab_filename(SPI_file);
end

if any(strcmp({'is','spinv','bin','mat'},ISformat)) & ~strcmp(SPIformat,'spi')
    if exist(IS_file,'file')
        SPI_file = '';
    else
        SPI_file = 'Solutionpoints (.spi)';
    end
end

if ~strcmp(SPIformat,'spi') & ~isempty(SPIformat)
    if ~isfield(settings,'method')
        settings.method = 'Gray';
    end
    if ~isfield(settings,'numspi')
        settings.numspi = 10000;
    end
    if ~isfield(settings,'atlasfile')
        settings.atlasfile = '';
    end
    if ~isfield(settings,'AtlasTemplate')
        settings.AtlasTemplate = '';
    end
    if ~isfield(settings,'RemoveBorder')
        settings.RemoveBorder = false;
    end
    if ~isfield(settings,'atlasexclude')
        settings.atlasexclude = [];
    end
    if ~isfield(settings,'ReduceAAL')
        settings.ReduceAAL = false;
    end
    
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Method for solutionpoints','method'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'Gray','GrayWhite','Surface','All'};
    Formats(end,1).span = [1 2];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Minimal number of solutionpoints','numspi'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'integer';
    Formats(end,1).limits = [0 99999]; % 4-digits (positive #)
    Formats(end,1).size = 50;
    Formats(end,1).span = [1 2];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Atlas-file','atlasfile'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.nii','MRI-file (.nii)'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).size = [300 0];
    
    Prompt(end+1,:) = {'Atlas Brain-Template','AtlasTemplate'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@get_template,'AtlasTemplate','AtlasTemplate','atlasfile'};
    
    Prompt(end+1,:) = {'Excluded regions','atlasexclude'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).limits = [0 9999]; % 4-digits (positive #)
    Formats(end,1).size = 350;
    Formats(end,1).callback = {@atlas_exclude,'@ALL','@ALL'};
    Formats(end,1).enable = 'inactive';
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Reduce AAL-Atlas to standard regions     ','ReduceAAL'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).format = 'integer';
    Formats(end,1).size = [-1 -1];
    Formats(end,1).callback = {@atlas_reduceaal,'@ALL','@ALL'};
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'Remove borders of regions     ','RemoveBorder'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).size = [-1 -1];
    Formats(end,1).span = [1 2];
    
    [settings,Cancelled] = inputsdlg(Prompt,'Create SPI',Formats,settings);
    if Cancelled == 1
        settings = [];
    end
else
    settings = [];
end

end

function settings = atlas_exclude(settings)
    if isempty(settings.atlasfile) | ~exist(settings.atlasfile,'file')
        return
    end
    atlas = lab_read_mri(settings.atlasfile);
    tmp = setdiff(unique(atlas.anatomy),0);
    if length(tmp) == 116
        list = lab_get_aal;
    else
        for i = 1:length(tmp)
            list{i,1} = ['Region ' num2str(i)];
        end
    end
    list = cat(1,cellstr('-none-'),list);
    if isempty(settings.atlasexclude)
        selection = 1;
    else
        selection = settings.atlasexclude + 1;
    end
    selection = listdlg('Name','Excluded regions','ListString',list,'SelectionMode','multiple', ...
        'ListSize',[200 400],'InitialValue',selection);
    selection = setdiff(selection,1);
    if ~isempty(selection)
        settings.atlasexclude = selection - 1;
    else
        settings.atlasexclude = [];
    end
end

function settings = atlas_reduceaal(settings)

if settings.ReduceAAL == true
    settings.ReduceAAL = false;
    settings.atlasexclude = [];
elseif ~isempty(settings.atlasfile) & exist(settings.atlasfile,'file')
    atlas = load_untouch_nii(settings.atlasfile);
    tmp = setdiff(unique(atlas.img),0);
    if length(tmp) == 116
        settings.ReduceAAL = true;
        settings.atlasexclude = setdiff(1:length(tmp),lab_get_AALselection);
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