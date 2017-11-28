function [settings] = lab_get_isfile(settings,isfiles,doshort)

SPI_file_Backup = settings.SPI_file;

if ~exist('settings','var') | ~isfield(settings,'IS_file') | strcmp(settings.IS_file,'Select File')
    if exist('isfiles','var') & isfield(isfiles,'IS_file')
        settings.IS_file = isfiles.IS_file;
    end
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt{end+1,1} = 'Inverse solution file';
    Formats(end+1,1).type = 'text';

    Prompt(end+1,:) = {'','IS_file'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.is','Cartool IS (.is)';'*.spinv','sLoreta IS (.spinv)'; ...
        '*.bin','Leadfiled (.bin)';'*.mat','Headmodel (.mat)'; ...
        '*.hdr;*.nii;*.fif;*.fiff','MRI-file (.hdr/.nii/.fif)'; ...
        '*.img;*.ima;*.dcm','DICOM (.img/.ima/.dcm)'};
    Formats(end,1).limits = [0 1]; % single file get
    Formats(end,1).size = [300 0];
        
    [settings,Cancelled] = inputsdlg(Prompt,'Inverse solution file',Formats,settings);
    if isempty(settings) | Cancelled == 1
        settings.IS_file = '';
        return
    elseif ~isempty(settings.IS_file) & exist(settings.IS_file,'file')
        isfiles = [];
        settings.IndividualFiles = false;
        isfiles.IS_file = settings.IS_file;
        isfiles = lab_find_isfiles([],[],isfiles);
    else
        settings.IndividualFiles = true;
    end
elseif ~isempty(settings.IS_file) & ~exist(settings.IS_file,'file')
    settings.IndividualFiles = true;
end


settings2 = lab_ISconvertnames2(settings);
[~,~,ISformat] = lab_filename(settings2.IS_file);
if isempty(ISformat)
    ISformat = settings.IS_file;
end
clearvars settings2

if settings.IndividualFiles == false & isfield(isfiles,'SPI_file')
    settings.SPI_file = isfiles.SPI_file;
    settings.LOC_file = isfiles.LOC_file;
    settings.MRI_file = isfiles.MRI_file;
    settings.ROIS_file = isfiles.ROIS_file;
end

switch ISformat
    case 'is'
        settings.type = 'Cartool';
        if ~isfield(settings,'issettings') | ~isfield(settings.issettings,'regularization') | ...
                isempty(settings.issettings.regularization)
            settings.issettings.regularization = 6;
        end
        if settings.IndividualFiles == true
            settings.SPI_file = 'Solutionpoints (.spi)';
            settings.LOC_file = 'Electrodes-file (.xyz)';
            settings.MRI_file = 'MRI-file (.hdr)';
            settings.ROIS_file = 'ROIS-file (.rois)';
        end
        if isfield(settings,'excludebad')
            settings.excludebad = false;
        end
    case 'spinv'
        settings.type = 'sLoreta';
        settings.issettings = [];
        if settings.IndividualFiles == true
            settings.SPI_file = 'Solutionpoints (.spi)';
            if ~any(strcmp({'Electrodes-file (.xyz)','Electrodes-file (.els)', ...
                    'Locations in input file'},settings.LOC_file))
                settings.LOC_file = 'Electrodes-file (.xyz)';
            end
            settings.MRI_file = 'MRI-file (.hdr)';
            settings.ROIS_file = 'ROIS-file (.rois)';
        end
        if isfield(settings,'excludebad')
            settings.excludebad = false;
        end
    case 'bin'
        if ~isfield(settings,'type') | isempty(settings.type) | strcmp(settings.type,'Cartool') | strcmp(settings.type,'sLoreta')
            settings.type = 'LCMV Beamformer';
            settings.issettings = [];
        end
        if settings.IndividualFiles == true
            settings.SPI_file = 'Solutionpoints (.spi)';
            if ~any(strcmp({'Electrodes-file (.xyz)','Electrodes-file (.els)', ...
                    'Locations in input file'},settings.LOC_file))
                settings.LOC_file = 'Electrodes-file (.xyz)';
            end
            settings.MRI_file = 'MRI-file (.hdr)';
            settings.ROIS_file = 'ROIS-file (.rois)';
        end
    case 'mat'
        if ~isfield(settings,'type') | isempty(settings.type) | strcmp(settings.type,'Cartool') | strcmp(settings.type,'sLoreta')
            settings.type = 'LCMV Beamformer';
            settings.issettings = [];
        end
        if settings.IndividualFiles == true
            settings.SPI_file = 'Solutionpoints (.spi)';
            if ~any(strcmp({'Electrodes-file (.xyz)','Electrodes-file (.els)', ...
                    'Locations in input file'},settings.LOC_file))
                settings.LOC_file = 'Electrodes-file (.xyz)';
            end
            settings.MRI_file = 'MRI-file (.hdr)';
            settings.ROIS_file = 'ROIS-file (.rois)';
        end
    otherwise
        if ~isfield(settings,'type') | isempty(settings.type) | strcmp(settings.type,'Cartool') | strcmp(settings.type,'sLoreta')
            settings.type = 'LCMV Beamformer';
            settings.issettings = [];
        end
        if settings.IndividualFiles == true
            if ~any(strcmp({'MRI-file (IS-file)','MRI-file (spi.hdr)', ...
                    'MRI-file (spi.nii)'},settings.SPI_file))
                settings.SPI_file = 'MRI-file (IS-file)';
                settings.SPI.method = 'Gray';
                settings.SPI.numspi = 10000;
            end
            if strcmp(ISformat,'fif') | strcmp(ISformat,'fiff')
                settings.LOC_file = 'Locations in input file';
            elseif ~any(strcmp({'Electrodes-file (.xyz)','Electrodes-file (.els)', ...
                    'Locations in input file'},settings.LOC_file))
                settings.LOC_file = 'Electrodes-file (.xyz)';
            end
            if ~any(strcmp({'MRI-Atlas (.nii)','MRI-Atlas (.hdr)'},settings.ROIS_file))
                settings.ROIS_file = 'MRI-Atlas (.nii)';
            end
            settings.MRI_file = 'MRI-file (IS-file)';
        end
end

if ~any(strcmp({'is','spinv','bin','mat'},ISformat))
    if ~isfield(settings,'LF')
        settings.LF = [];
    end
    if ~isfield(settings,'MRI')
        settings.MRI = [];
    end
    settings.MRI = lab_get_pMRI(settings.MRI,settings.IS_file,settings.LF,isfiles);
    settings.LF.MRI = settings.MRI;
    settings.LF = lab_get_LF(settings.LF,ISformat,settings.LOC_file);
else
    settings.LF = [];
    settings.MRI = [];
end

if (isfield(settings,'resultscalar') & settings.resultscalar == 1) | ~any(strcmp({'is','spinv'},ISformat))
    if isfield(settings,'COV') & isempty(settings.COV)
        settings.COV = lab_get_cov(settings.COV);
    end
elseif isfield(settings,'COV')
    settings.COV = [];
end

[~,~,tmp] = lab_filename(settings.SPI_file);
if strcmp(tmp,'hdr') | strcmp(tmp,'nii')
    [settings.SPI,settings.SPI_file] = lab_get_spi(settings.SPI,settings.SPI_file,settings.IS_file);
end

if strcmp(SPI_file_Backup,'None')
    settings.SPI_file = 'None';
    settings.SPI = [];
end