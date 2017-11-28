function settings = lab_get_issettings(settings)

    
settings2 = lab_ISconvertnames2(settings);
[~,~,ISformat] = lab_filename(settings2.IS_file);
if isempty(ISformat)
    ISformat = settings.IS_file;
end
switch ISformat
    case 'is'
        settings.type = 'Cartool';
    case 'spinv'
        settings.type = 'sLoreta';
end
clearvars settings2 ISformat

Prompt = cell(0,2);
Formats = [];

switch settings.type
    case 'Cartool'
        if ~isfield(settings,'issettings') | ~isfield(settings.issettings,'regularization') | ...
                isempty(settings.issettings.regularization)
            settings.issettings.regularization = 6;
        end
        Prompt(end+1,:) = {'Regularization','regularization'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).limits = [0 12];
        Formats(end,1).size = 25;
    case 'sLoreta'
        settings.issettings = [];
    case 'eLoreta'
        settings.issettings = [];
    case 'LCMV Beamformer'
        if ~isfield(settings,'issettings') | ~isfield(settings.issettings,'wn')
            settings.issettings.wn = false;
            settings.issettings.LCMVcn = true;
            settings.issettings.progressbar = false;
        end
        Prompt(end+1,:) = {'Normalize Leadfield','LCMVcn'};
        Formats(end+1,1).type = 'check';
        
        Prompt(end+1,:) = {'Normalize Weights','wn'};
        Formats(end+1,1).type = 'check';
        
        Prompt(end+1,:) = {'Show Progressbar','progressbar'};
        Formats(end+1,1).type = 'check';
    case 'Minimum Norm Estimate'
        if ~isfield(settings,'issettings') | ~isfield(settings.issettings,'lambda') | ...
                isempty(settings.issettings.lambda)
            settings.issettings.lambda = 3;
            settings.issettings.prewhiten = true;
            settings.issettings.scalesourcecov = true;
        end
        Prompt(end+1,:) = {'Lambda','lambda'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).limits = [0 20];
        Formats(end,1).size = 25;
        
        Prompt(end+1,:) = {'Prewhiten leadfield','prewhiten'};
        Formats(end+1,1).type = 'check';

        Prompt(end+1,:) = {'Scale source covariance matrix','scalesourcecov'};
        Formats(end+1,1).type = 'check';
    case 'Dipol fit'
        if ~isfield(settings,'issettings') | ~isfield(settings.issettings,'numdipols') | ...
                isempty(settings.issettings.numdipols)
            settings.issettings.numdipols = 1;
            settings.issettings.model = 'regional';
            settings.issettings.symmetry = '';
        end
        settings.issettings.SPI_file = settings.SPI_file;
        settings.issettings.SPI = settings.SPI;
        
        Prompt(end+1,:) = {'Number of dipoles','numdipols'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).limits = [0 100];
        Formats(end,1).size = 25;
        Formats(end,1).callback = {@set_numdipols,'@ALL','@ALL'};
        Formats(end,1).span = [1 2];
        
        Prompt(end+1,:) = {'Model','model'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).format = 'input';
        Formats(end,1).items = {'regional','moving'};
        Formats(end,1).span = [1 2];
        
        Prompt(end+1,:) = {'Symmetry (number = 2)','symmetry'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).format = 'input';
        Formats(end,1).items = {'','x','y','z'};
        Formats(end,1).callback = {@set_symmetry,'@ALL','@ALL'};
        Formats(end,1).span = [1 2];
        
        Prompt(end+1,:) = {'Grid search','SPI_file'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).format = 'input';
        Formats(end,1).items = {'None','Select File','Solutionpoints (.spi)','MRI-file (IS-file)','MRI-file (spi.hdr)','MRI-file (spi.nii)'};
        Formats(end,1).size = 300;
        Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file','IS_file'};
        
        Prompt(end+1,:) = {'settings','SPI'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).callback = {@lab_get_spi,{'SPI','SPI_file'},'SPI','SPI_file',settings.IS_file};
end

if ~isempty(Formats)
    [settings.issettings,Cancelled] = inputsdlg(Prompt,'IS settings',Formats,settings.issettings);
    if Cancelled == 1
        settings.issettings = [];
    elseif isfield(settings.issettings,'SPI_file')
        if ~strcmp(settings.issettings.SPI_file,'None')
            settings.SPI_file = settings.issettings.SPI_file;
            settings.SPI = settings.issettings.SPI;
        else
            settings.SPI_file = '';
            settings.SPI = [];
        end
    end
end

end

function settings = set_numdipols(settings)
   if settings.numdipols > 1
       settings.gridsearch = false;
   end
   if settings.numdipols ~= 2
       settings.symmetry = '';
   end
end

function settings = set_symmetry(settings)
   if ~isempty(settings.symmetry)
       settings.gridsearch = false;
       settings.numdipols = 2;
   end
end