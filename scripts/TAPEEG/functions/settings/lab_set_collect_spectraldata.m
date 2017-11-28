function [cfg,Files,skipprocessing] = lab_set_collect_spectraldata(cfg,nofolder)

disp ('   Ask for collect spectraldata settings')
    
skipprocessing = 0;
if ~exist('nofolder','var')
    nofolder = false;
end
Files = {};
Patient = ' ';
Filename = ' ';

if ~exist('cfg','var') | ~isfield(cfg,'CollectFFT') | ~isfield(cfg.CollectFFT,'source')
    cfg.CollectFFT.searchfolder = '';
    cfg.CollectFFT.includestring{1} = '';
    cfg.CollectFFT.excludestring{1} = '';
    cfg.CollectFFT.searchIS = false;
    cfg.CollectFFT.subjectname = 0;
    cfg.CollectFFT.subjecttext = ' ';
    cfg.CollectFFT.folder = 'FrequencyAnalysis';
    cfg.CollectFFT.mappings = [];
    cfg.CollectFFT.Locs = [];
    cfg.CollectFFT.mappingBA = [];
    cfg.CollectFFT.lowfreqpeak = 4;
    cfg.CollectFFT.highfreqpeak = 14;
    cfg.CollectFFT.lowfreqcog = 4;
    cfg.CollectFFT.highfreqcog = 14;
    cfg.CollectFFT.source = 'median';
    cfg.CollectFFT.MinPeak2Min = 1.3;
    cfg.CollectFFT.spectralbands = lab_get_spectralbands;
    cfg.CollectFFT.spectralbandsI = false;
    cfg.CollectFFT.IndivBands = false;
    cfg.CollectFFT.correctpf = true;
    cfg.CollectFFT.calcsingle = true;
    cfg.CollectFFT.qualityrange = 0.5;
    cfg.CollectFFT.qualityplot = true;
    cfg.CollectFFT.QUALITY = [];
    cfg.CollectFFT.plotchans = false;
    cfg.CollectFFT.plotmappings = false;
end
if isempty(cfg.CollectFFT.mappings)
    if (~isfield(cfg.CollectFFT,'searchIS') | cfg.CollectFFT.searchIS == false) & ...
            isfield(cfg,'FFT') & isfield(cfg.FFT,'mappings')
        cfg.CollectFFT.mappings = cfg.FFT.mappings;
    elseif isfield(cfg,'ISFFT') & isfield(cfg.ISFFT,'mappings')
         cfg.CollectFFT.mappings = cfg.ISFFT.mappings;
    end
end
if nofolder == true
    cfg.CollectFFT.searchfolder = '';
    cfg.CollectFFT.Files = {};
elseif ~isfield(cfg.CollectFFT,'Files')
    cfg.CollectFFT.Files = {};
end

doQUALITYchans = true;
doQUALITYact = true;
doQUALITYepoch = true;
Channels = [];

Prompt = cell(0,2);
Formats = [];

if nofolder == false
    Prompt{end+1,1} = 'Search folder';
    Formats(end+1,1).type = 'text';
    Formats(end,1).span = [1 4];
    
    Prompt(end+1,:) = {'','searchfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'dir';
    Formats(end,1).size = 300;
    Formats(end,1).span = [1 3];
    Formats(end,1).callback = {@set_folder,'@ALL','@ALL'};
end

Prompt{end+1,1} = 'Include strings';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'','includestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'ISFFT-Results','searchIS'};
Formats(end+1,1).type = 'check';

Prompt{end+1,1} = 'Exclude strings';
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'','excludestring'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 300;
Formats(end,1).span = [1 4];

if nofolder == false
    Prompt(end+1,:) = {'Search Files','Files'};
    Formats(end+1,1).type = 'check';
    Formats(end,1).callback = {@search_files,'@ALL','@ALL','$subjecttext','$patienttext'};
end

Prompt(end+1,:) = {'Preload File-Info',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [110 25];
Formats(end,1).callback = {@load_fileinfo,'@ALL','@ALL','$subjecttext','$patienttext'};
Formats(end,1).span = [1 3];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

if isfield(cfg.CollectFFT,'subjecttext');
    Prompt(end+1,:) = {cfg.CollectFFT.subjecttext,'subjecttext'};
else
    Prompt(end+1,:) = {'','subjecttext'};
end
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Number of underscores in subject name','subjectname'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [-99 99];
Formats(end,1).size = 30;
Formats(end,1).callback = {@set_patient,'@ALL','@ALL','$patienttext'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {Patient,'patienttext'};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 2];

Formats(end+1,1).type = 'none';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Background Channels','mappingBA'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).enable = 'inactive';
Formats(end,1).callback = {@load_background,'mappingBA','mappingBA','Locs'};
Formats(end,1).size = 350;
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Output-folder', 'folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 150;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'LOCS','Locs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_locs,'Locs','Locs'};

Prompt(end+1,:) = {'Mappings','mappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_mappings,'mappings','mappings','Locs'};

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Source for peak and median freq','source'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'mean','median'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Spectral Bands','spectralbands'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = 60;
Formats(end,1).callback = {@lab_get_spectralbands,'spectralbands','spectralbands','spectralbandsI',true};

Prompt(end+1,:) = {'Individual Bands','spectralbandsI'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'input';
Formats(end,1).callback = {@do_spectralbandsI,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Range peak frequency (Hz)    low:','lowfreqpeak'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'high:','highfreqpeak'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Range median frequency (Hz) low:','lowfreqcog'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'high:','highfreqcog'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Minimal Peak2Min','MinPeak2Min'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_MinPeak2Min,'MinPeak2Min','MinPeak2Min'};
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Calculate single channel peaks','calcsingle'};
Formats(end+1,1).type = 'check';
Formats(end,1).format = 'integer';
Formats(end,1).span = [1 4];
Formats(end,1).callback = {@check_single,{'calcsingle','calcL2Rratio'},'calcsingle','calcL2Rratio'};

Prompt(end+1,:) = {'Manually correct peak frequency','correctpf'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Quality','qualityplot'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Range (+/-Hz)','qualityrange'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Evaluate','QUALITY'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_quality,'QUALITY','QUALITY'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Plot spectras for channels','plotchans'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Plot spectras for mappings','plotmappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 2];

[cfg.CollectFFT,Cancelled] = inputsdlg(Prompt,'Collect spectral data',Formats,cfg.CollectFFT);
if isempty(cfg.CollectFFT) | Cancelled == 1
    cfg.CollectFFT = [];
    skipprocessing = 1;
    return
else
    pause(0.2);
    Files = cfg.CollectFFT.Files;
    cfg.CollectFFT = rmfield(cfg.CollectFFT,'Files');
    if isfield(cfg.CollectFFT,'searchfolder') & exist(cfg.CollectFFT.searchfolder,'dir')
        save(fullfile(cfg.CollectFFT.searchfolder,'settings_spectras.mat'),'cfg','-v7.3');
    end
end

    function settings = set_folder(settings)
        if exist(fullfile(settings.searchfolder,'settings_spectras.mat'),'file')
            searchfolder = settings.searchfolder;
            TMP = load(fullfile(settings.searchfolder,'settings_spectras.mat'));
            if isfield(TMP,'cfg') & isfield(TMP.cfg,'CollectFFT') & ~isempty(TMP.cfg.CollectFFT)
                settings = TMP.cfg.CollectFFT;
                if isfield(settings,'spectralbands') & isnumeric(settings.spectralbands) & size(settings.spectralbands,2) == 2
                    settings.spectralbands = cat(2,cell(size(settings.spectralbands,1),1),num2cell(settings.spectralbands));
                end
                if ~isfield(settings,'spectralbandsI')
                    settings.spectralbandsI = false;
                end
                settings.searchfolder = searchfolder;
            end
            clearvars TMP searchfolder
        end
    end
    
    function settings = load_fileinfo(settings,Shandle,Phandle)
        TFiles = settings.Files;
        if isempty(TFiles)
            [Filename,Filepath] = uigetfile('*.mat','Select *Spect.mat');
            TFiles{1} = fullfile(Filepath,Filename);
        end
        Subjectname = {'',''};
        if ~isempty(TFiles);
            Filename = TFiles{1};
            try %#ok<TRYNC>
                MAT = load(Filename);
                Subjectname = lab_prepare_subjectname(Filename);
                if isfield(MAT,'header') & isfield(MAT.header,'patient')
                    Patient = MAT.header.patient;
                    set(Phandle,'String',['= ' regexprep(Patient,'_',' ')]);
                    if isempty(settings.subjectname)
                        settings.subjectname = [];
                    end
                end
                if isfield(MAT,'header')
                    Channels = cellstr(MAT.header.channels);
                    if isfield(MAT.header,'numdatachannels')
                        Channels = Channels(1:MAT.header.numdatachannels,1);
                    end
                    if isfield(MAT.header,'locs')
                        settings.Locs = MAT.header.locs;
                    end
                    if isfield(MAT.header,'badchans')
                        doQUALITYchans = true;
                    end
                    if isfield(MAT.header,'activationsexcluded')
                        doQUALITYact = true;
                    end
                    if isfield(MAT.header,'quality')
                        doQUALITYepoch = true;
                    end
                    clearvars MAT
                end
            end
            if isempty(settings.mappingBA) & ~isempty(Channels)
                settings.mappingBA = lab_get_bactivity(size(Channels,1));
            end
        else
            settings.mappingBA = [];
            settings.folder = [];
            settings.Locs = [];
            settings.mappings = [];
        end
        if ~isempty(Subjectname{1})
            settings.subjecttext = [Subjectname{1} ' ' Subjectname{2}];
        else
            settings.subjecttext = Subjectname{2};
        end
        set(Shandle,'String',settings.subjecttext);
    end
    
    function settings = search_files(settings,Shandle,Phandle)
        if isempty(settings.Files)
            cfg2.CollectFFT = settings;
            TFiles = lab_collect_spectraldata_search(cfg2);
            if ~isempty(TFiles)
                settings.Files = TFiles;
            end
        else
            disp ('Select Files')
            selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
                'ListString',settings.Files,'InitialValue',1:length(settings.Files),'CancelString','None','ListSize',[450 400]);
            pause(0.2);
            if ~isempty(selection)
                settings.Files = settings.Files(1,selection);
            else
                settings.Files = {};
            end
        end
        if ~isempty(settings.Files)
            settings = load_fileinfo(settings,Shandle,Phandle);
        end
    end
    
    function settings = set_patient(settings,Phandle)
        if ~isempty(settings.subjectname)
            set(Phandle,'String',['= ' regexprep(lab_subjectname(Filename,settings),'_',' ')]);
        else
            set(Phandle,'String',['= ' regexprep(Patient,'_',' ')]);
        end
    end
    
    function Mappings = load_mappings(Mappings,LOCS)
        settings2 = cfg;
        if ~isempty(Channels)
            settings2.EXTRA.numdatachans = size(Channels,1);
        end
        Mappings = lab_load_mappings(Mappings,settings2,[],LOCS);
    end
    
    function Background = load_background(Background,LOCS)
        Background = lab_load_background(Background,Channels,LOCS);
    end
    
    function LOCS = load_locs(LOCS)
        LOCS = lab_load_locs(LOCS,cfg,size(Channels,1));
    end
    
    function QUALITY = get_quality(QUALITY)
        QUALITY = lab_get_QUALITY(QUALITY,0,doQUALITYchans,doQUALITYact,doQUALITYepoch);
    end
    
    function [calcsingle,calcL2Rratio] = check_single(calcsingle,calcL2Rratio)
        if calcsingle == false
            calcsingle = true;
        else
            calcsingle = false;
            calcL2Rratio = false;
        end
    end
    
    function MinPeak2Min = set_MinPeak2Min(MinPeak2Min)
        if isempty(MinPeak2Min)
            MinPeak2Min = 1.3;
        end
        settings.MinPeak2Min = MinPeak2Min;
        
        Prompt2 = {'Minimal Peak2Min for detection of peak frequency','MinPeak2Min'};
        Formats2.type = 'edit';
        Formats2.format = 'float';
        Formats2.limits = [0 inf];
        Formats2.size = 30;
        
        [settings,Cancelled2] = inputsdlg(Prompt2,'Minimal Peak2Min',Formats2,settings);
        if isempty(settings) | Cancelled2 == 1
            MinPeak2Min = 1.3;
        else
            MinPeak2Min = settings.MinPeak2Min;
        end
        clearvars Prompt2 Formats2 Cancelled2 settings
    end
end

function settings = do_spectralbandsI(settings)
    if settings.spectralbandsI == false
        settings.spectralbandsI = true;
        settings.spectralbands = [1,2,4,5,6,8,9];
    else
        settings.spectralbandsI = false;
        settings.spectralbands = lab_get_spectralbands;
    end
end



