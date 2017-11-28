function ROIS = lab_load_ROIS(ROIS,SPI_file,SPI)

if ~exist('SPI','var')
    SPI = [];
end
if ~exist('SPI_file','var')
    SPI_file = [];
    SPIformat = '';
elseif ischar(SPI_file)
    [~,~,SPIformat] = lab_filename(SPI_file);
else
    SPIformat = '';
end
if ~exist('ROIS','var')
    settings.rois = [];
else
    settings.rois = ROIS;
end
settings.ROIS_file = '';
settings.ROIS = [];

if strcmp(SPIformat,'hdr') | strcmp(SPIformat,'nii')
    Filetype = {'*.hdr;*.nii','MRI-file (.hdr/.nii)'};
elseif strcmp(SPIformat,'spi')
    Filetype = {'*.rois','ROIS-file (.rois)'};
else
    Filetype = {'*.rois','ROIS-file (.rois)';'*.hdr;*.nii','MRI-file (.hdr/.nii)'};
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'ROIS - File','ROIS_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = Filetype;
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).callback =  {@read_rois,{'ROIS','ROIS_file'},'ROIS_file','ROISsettings',SPI_file,SPI};
Formats(end,1).size = 270;

Prompt(end+1,:) = {'Settings','ROISsettings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_ROIS,'@ALL','@ALL'};

Prompt(end+1,:) = {'ROIS','ROIS'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [190 90];
Formats(end,1).enable = 'inactive';

[settings,Cancelled] = inputsdlg(Prompt,'Inverse solution',Formats,settings);
if isempty(settings) | Cancelled == 1
    ROIS = [];
else
    ROIS = settings.ROIS;
end

end

function [ROIS,ROIS_file] = read_rois(ROIS_file,ROISsettings,SPI_file,SPI)
   if strcmp(SPI_file(end-3:end),'.spi')
       slocs = lab_read_spi(SPI_file);
   else
       slocs = lab_create_sp(SPI_file,SPI);
   end
   ROIS = lab_create_rois(ROIS_file,ROISsettings,slocs,SPI_file);
   if isempty(ROIS)
       ROIS_file = '';
   end
end
