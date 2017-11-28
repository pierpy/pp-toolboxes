function [settings,Prompt,Formats] = lab_get_isfiles(settings,doshort,doisfft,isfiles,Prompt,Formats)

if ~exist('isfiles','var')
    isfiles = [];
end
if ~exist('doshort','var')
    doshort = 0;
end
if ~exist('doisfft','var')
    doisfft = 0;
end
if ~exist('Formats','var')
    if exist('settings','var') & ~isempty(settings) & isfield(settings,'IS_file')
        settings = lab_ISconvertnames(settings);
    else
        settings = [];
    end
    Prompt = cell(0,2);
    Formats = [];
    dodiag = 1;
else
    dodiag = 0;
end
if ~isfield(settings,'IS_file')
    settings.IS_file = ' ';
    settings.SPI_file = ' ';
    settings.SPI =[];
    settings.LOC_file = ' ';
    settings.MRI_file = ' ';
    settings.ROIS_file = ' ';
    settings.ROISmethods = {};
    settings.spectralbands = [];
    settings.spectralbandsI = false;
    settings.ROISsettings = [];
    settings.LF = [];
    settings.COV = [];
    settings.IndividualFiles = false;
end

Prompt(end+1,:) = {'Input-file for Inverse solution','IS_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Select File','Cartool (.is)','sLoreta (.spinv)','Leadfield (LF.bin)', ...
    'Headmodel (.mat)','MRI-file (.hdr)','MRI-file (.nii)','MRI-file (iso.fif)', ...
    'MRI-file (mri.fif)','MRI-file (dicom)'};
Formats(end,1).size = 300;
Formats(end,1).callback = {@lab_get_isfile,'@ALL','@ALL',isfiles,doshort};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Individual Files','IndividualFiles'};
Formats(end+1,1).type = 'check';
Formats(end,1).enable = 'inactive';

Prompt(end+1,:) = {'MRI settings','MRI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_pMRI,{'MRI','LF'},'MRI','IS_file','LF',isfiles};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Leadfield','LF'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_LF,'LF','LF','IS_file','LOC_file'};

if isfield(settings,'resultscalar')
    Prompt(end+1,:) = {'Covariance matrix','COV'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_cov,'COV','COV'};
    Formats(end,1).span = [1 2];
else
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 2];
end

Prompt(end+1,:) = {'Type of inverse solution','type'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
if doisfft == 0
    Formats(end,1).items = {'','LCMV Beamformer','sLoreta','eLoreta','Minimum Norm Estimate','Dipol fit'};
else
    Formats(end,1).items = {'','LCMV Beamformer','sLoreta','eLoreta','Minimum Norm Estimate'};
end
Formats(end,1).callback = {@lab_get_issettings,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Settings','issettings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_issettings,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Type of Solutionpoints file','SPI_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Select File','Solutionpoints (.spi)','MRI-file (IS-file)','MRI-file (spi.hdr)','MRI-file (spi.nii)'};
Formats(end,1).size = 300;
Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file','IS_file'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'create SPI','SPI'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file','IS_file'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Type of Electrodes file','LOC_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'Select File','Locations in input file','Electrodes-file (.xyz)','Electrodes-file (.els)'};
Formats(end,1).size = 300;
Formats(end,1).callback = {@lab_get_LOC,'LOC_file','LOC_file'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 3];

if doshort == 0
    Prompt(end+1,:) = {'MRI-file (visualisation)','MRI_file'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'Select File','MRI-file (IS-file)','MRI-file (.hdr)'};
    Formats(end,1).size = 300;
    Formats(end,1).callback = {@lab_get_MRI,'MRI_file','MRI_file'};
    Formats(end,1).span = [1 3];
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
    
    Prompt(end+1,:) = {'ROIS-file','ROIS_file'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'Select File','MRI-Atlas (.nii)','MRI-Atlas (.hdr)','ROIS-file (.rois)','LOC-file (.xyz)'};
    Formats(end,1).size = 220;
    Formats(end,1).callback = {@lab_get_ROIS,'@ALL','@ALL'};
    Formats(end,1).span = [1 2];
    
    Prompt(end+1,:) = {'settings','ROISsettings'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@lab_get_ROIS,'@ALL','@ALL'};
    
    Prompt(end+1,:) = {'ROIS method','ROISmethods'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'listbox';
    Formats(end,1).format = 'input';
    if doisfft == 1
        Formats(end,1).items = {'mean';'maxpower';'center';'pseudoelectrodes';'off'};%G Bogaarts edit 2016
        Formats(end,1).size = [110 80];
    else
        Formats(end,1).items = {'maxpower';'center';'pseudoelectrodes';'off'};%G Bogaarts edit 2016
        Formats(end,1).size = [110 65];
        Formats(end,1).callback = {@set_roismethods,'@ALL','@ALL'};
    end
    Formats(end,1).limits = [0 2]; % multi-select
    Formats(end,1).span = [2 1];
    
    if doisfft == 0
        Prompt(end+1,:) = {'Spectral Bands','spectralbands'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).callback = {@do_spectralbands,'spectralbands','spectralbands','spectralbandsI','ROISmethods'};
        Formats(end,1).span = [1 2];
        
        Prompt(end+1,:) = {'Individual Bands','spectralbandsI'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).format = 'input';
        Formats(end,1).callback = {@do_spectralbandsI,{'spectralbands','spectralbandsI'},'spectralbands','spectralbandsI','ROISmethods'};
        Formats(end,1).span = [1 2];
    else
        Formats(end+1,1).type = 'none';
        Formats(end,1).span = [2 2];
    end
    
    Formats(end+1,1).type = 'none';
    Formats(end,1).span = [1 3];
end

if dodiag == 1
    [settings,Cancelled] = inputsdlg(Prompt,'Inverse solution',Formats,settings);
    if isempty(settings) | Cancelled == 1
        settings = [];
    else
        settings = lab_ISconvertnames2(settings);
    end
end

end

function settings = set_roismethods(settings)
   if any(strcmp(settings.ROISmethods,'pseudoelectrodes'))
       settings.spectralbands = [];
       settings.spectralbandsI = false;
       if ~isempty(settings.ROIS_file)
           if strcmp(settings.ROIS_file,'MRI-Atlas (.nii)') | strcmp(settings.ROIS_file,'MRI-Atlas (.hdr)') | ...
                   strcmp(settings.ROIS_file,'ROIS-file (.rois)')
               settings.ROIS_file = 'LOC-file (.xyz)';
           elseif ~strcmp(settings.ROIS_file,'LOC-file (.xyz)')
               [~,~,ROIS_format] = lab_filename(settings.ROIS_file);
               if ~strcmp(ROIS_format,'xyz') & ~strcmp(ROIS_format,'els')
                   settings.ROIS_file = 'Select File';
                   settings = lab_get_ROIS(settings);
               end
           end
       end
   else
       if any(strcmp(settings.ROISmethods,'maxpower'))
           if isempty(settings.spectralbands)
               settings.spectralbands = lab_get_spectralbands(settings.spectralbands,settings.spectralbandsI);
           end
       else
           settings.spectralbands = [];
       end
       if ~isempty(settings.ROIS_file)
           if strcmp(settings.ROIS_file,'LOC-file (.xyz)')
               settings.ROIS_file = '';
           elseif ~strcmp(settings.ROIS_file,'MRI-Atlas (.nii)') & ~strcmp(settings.ROIS_file,'MRI-Atlas (.hdr)') & ...
                   ~strcmp(settings.ROIS_file,'ROIS-file (.rois)')
               [~,~,ROIS_format] = lab_filename(settings.ROIS_file);
               if ~strcmp(ROIS_format,'nii') & ~strcmp(ROIS_format,'hdr') & ~strcmp(ROIS_format,'rois')
                   settings.ROIS_file = 'Select File';
                   settings = lab_get_ROIS(settings);
               end
           end
       end
   end
end

function spectralbands = do_spectralbands(spectralbands,spectralbandsI,ROISmethods)
    if any(strcmp(ROISmethods,'maxpower'))
        spectralbands = lab_get_spectralbands(spectralbands,spectralbandsI,true);
    else
        spectralbands = [];
    end
end

function [spectralbands,spectralbandsI] = do_spectralbandsI(spectralbands,spectralbandsI,ROISmethods)
    if ~any(strcmp(ROISmethods,'maxpower'))
        spectralbandsI = false;
        spectralbands = [];
    elseif spectralbandsI == false
        spectralbandsI = true;
        if ~isnumeric(spectralbands) | isempty(spectralbands)
            spectralbands = [1,2,4,5,6,8,9];
        end
    else
        spectralbandsI = false;
        if isnumeric(spectralbands) | isempty(spectralbands)
            spectralbands = lab_get_spectralbands;
        end
    end
end