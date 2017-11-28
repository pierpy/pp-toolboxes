function lab_detect_dbs(cfg)

disp('Detect DBS-stimulator artifacts')

if ~exist('cfg','var')
    cfg = [];
end

disp('Select mode')
button = questdlg('Single File or Folder','Select mode','Cancel','Folder','File','File');
if strcmp(button,'Folder')
    [calc,cfg] = lab_search_files;
elseif strcmp(button,'File')
    disp ('Ask for EEG-file')
    [EEG_file,EEG_filepath]=uigetfile('*.*','Select EEG file');
    if ~isempty(EEG_file) & EEG_file ~= 0
        cd(EEG_filepath);
        if exist(fullfile(EEG_filepath,'signal1.bin'),'file')
            tmp = strfind(EEG_filepath,filesep);
            EEG_file = [EEG_filepath(tmp(end-1)+1:tmp(end)-1) '.mff'];
            EEG_filepath = EEG_filepath(1:tmp(end-1));
            clearvars tmp
        elseif exist(fullfile(fullfile(EEG_filepath,EEG_file),'signal1.bin'),'file')
            EEG_file = [EEG_file '.mff'];
        end
        calc.Filelist{1,1} = fullfile(EEG_filepath,EEG_file);
        calc.Filelist_done = [];
        calc.Filelist_doneall = [];
        cfg.SEARCH.searchfolder = 0;
        clearvars EEG_file EEG_filepath
    else
        return
    end
else
    return
end

if ~isfield(cfg,'DBS') | isempty(cfg.DBS)
    [cfg,skipprocessing] = lab_set_detect_dbs(cfg,calc);
    if skipprocessing == 1
        return
    end
end

for L = 1:length(calc.Filelist);
    disp(['   Read ' calc.Filelist{L}])
    [data,header,cfg] = lab_read_data(calc.Filelist{L},cfg);
    if ~isempty(data)
        if ~isfield(header,'numdatachannels')
            header.numdatachannels = header.numchannels;
        end
        
        if ~isfield(cfg,'EEG_filepath')
            if isfield(header,'EEG_filepath')
                cfg.EEG_filepath = header.EEG_filepath;
            else
                cfg.EEG_filepath = pwd;
            end
        end
        [cfg.EEG_fileS,cfg.DBS] = lab_subjectname(calc.Filelist{L},cfg.DBS);
        cfg.EEG_file = [cfg.EEG_fileS '.sef'];
        
        % highpass filter
        settings.FILT.filtermode = 'butter';
        settings.FILT.filtorder = 2;
        settings.FILT.Filter.firstchan = 1;
        settings.FILT.Filter.lastchan = header.numdatachannels;
        settings.FILT.Filter.highpass = 1;
        settings.FILT.Filter.lowpass = [];
        settings.FILT.Filter.notch = [];
        [dataT,header] = lab_filtering(data,header,settings);
        
        % split in segments
        disp('   Find segments with high voltage artifact')
        Power = mean(dataT(1:header.numdatachannels,:),1).^2;
        clearvars dataT
        tmp = resample(Power,1,header.samplingrate);
        tmp = tmp > cfg.DBS.threshold;
        tmp2 = diff(tmp);
        
        % find high voltage segments
        StartS = find(tmp2 == 1) + 1;
        if tmp(1) == 1
            StartS = [1 StartS]; %#ok<AGROW>
        end
        StartS = StartS * header.samplingrate - (header.samplingrate-1);
        StopS = find(tmp2 == -1);
        if tmp(end) == 1
            StopS = [StopS length(tmp)]; %#ok<AGROW>
        end
        StopS = StopS * header.samplingrate;
        StopS(StopS > size(data,2)) = size(data,2);
        if length(StartS) == length(StopS)
            Seg = [StartS(:) StopS(:)];
            Seg = cat(2,Seg,ones(size(Seg,1),1));
        else
            Seg = [];
            disp('     no segments with high power artifact found')
        end
        
        % find low voltage segments
        StartS = find(tmp2 == -1) + 1;
        if tmp(1) == 0
            StartS = [1 StartS]; %#ok<AGROW>
        end
        StartS = StartS * header.samplingrate - (header.samplingrate-1);
        StopS = find(tmp2 == 1);
        if tmp(end) == 0
            StopS = [StopS length(tmp)]; %#ok<AGROW>
        end
        StopS = StopS * header.samplingrate;
        StopS(StopS > size(data,2)) = size(data,2);
        if length(StartS) == length(StopS)
            Seg2 = [StartS(:) StopS(:)];
            Seg2 = cat(2,Seg2,zeros(size(Seg2,1),1));
            Seg = cat(1,Seg,Seg2);
            clearvars Seg2
        end
        clearvars tmp tmp2 StartS StopS Power
        
        % sort and correct segments
        if size(Seg,1) > 1
            [~,Idx] = sort(Seg(:,1));
            Seg = Seg(Idx,:);
            Seg(:,4) = Seg(:,2) - Seg(:,1);
            tmp = find(Seg(:,4) <= 5 * header.samplingrate);
            Seg(tmp,3) = 0;
            for i = tmp(:)'
                if Seg(i,3) == 0
                    if i == 1
                        if Seg(i+1,3) == 1
                            Seg(i,3) = 1;
                        end
                    elseif i == size(Seg,1)
                        if Seg(i-1,3) == 1
                            Seg(i,3) = 1;
                        end
                    elseif Seg(i+1,3) == 1 | Seg(i-1,3) == 1
                        Seg(i,3) = 1;
                    end
                end
            end
            tmp = 1;
            while tmp < size(Seg,1)
                if Seg(tmp,3) == Seg(tmp+1,3)
                    Seg(tmp+1,1) = Seg(tmp,1);
                    Seg = Seg(setdiff(1:size(Seg,1),tmp),:);
                else
                    tmp = tmp + 1;
                end
            end
            clearvars i tmp Idx
            Seg(:,4) = Seg(:,2) - Seg(:,1);
        end
        
        if ~isfield(cfg.DBS,'folder')
            cfg.DBS.folder = 'RemoveDBS';
        end
        cfg.DBS.EEG_filepath = fullfile(cfg.EEG_filepath,cfg.DBS.folder);
        warning off %#ok<WNOFF>
        mkdir(cfg.DBS.EEG_filepath);
        warning on %#ok<WNON>
        
        % Store Stimulator on/off as marker-file
        if ~isfield(header,'events') | ~isfield(header.events,'POS')
            header.events.POS = [];
            header.events.OFF = [];
            header.events.DUR = [];
            header.events.TYP = {};
        end
        for i = 1:size(Seg,1)
            if Seg(i,3) == 1
                header.events.POS(end+1) = int64(Seg(i,1));
                header.events.OFF(end+1) = int64(0);
                header.events.DUR(end+1) = int64(1);
                header.events.TYP{end+1} = 'StimulatorOn';
                header.events.POS(end+1) = int64(Seg(i,2));
                header.events.OFF(end+1) = int64(0);
                header.events.DUR(end+1) = int64(1);
                header.events.TYP{end+1} = 'StimulatorOff';
            end
        end
        if ~isempty(header.events.POS)
            [header.events.POS,Idx] = sort(header.events.POS);
            header.events.DUR = header.events.DUR(1,Idx);
            header.events.OFF = header.events.OFF(1,Idx);
            header.events.TYP = header.events.TYP(1,Idx);
            lab_write_mrk(fullfile(cfg.DBS.EEG_filepath,[cfg.EEG_file '.mrk']),header)
        end
        
        DBS = 0;
        NoDBS = 0;
        for i = 1:size(Seg,1)
            headerS = header;
            headerS.numtimeframes = Seg(i,2) - Seg(i,1) + 1;
            if isfield(headerS,'events') & isfield(headerS.events,'POS') & ~isempty(headerS.events.POS)
                headerS.events.POS = headerS.events.POS - Seg(i,1) + 1;
                tmp = find(headerS.events.POS <= headerS.numtimeframes);
                tmp = intersect(tmp,find(headerS.events.POS > 0));
                headerS.events.POS = headerS.events.POS(1,tmp);
                headerS.events.DUR = headerS.events.DUR(1,tmp);
                headerS.events.OFF = headerS.events.OFF(1,tmp);
                headerS.events.TYP = headerS.events.TYP(1,tmp);
                clearvars tmp
            end
            if Seg(i,3) == 1
                header.dbs = true;
                DBS = DBS + 1;
                lab_write_sef(fullfile(cfg.DBS.EEG_filepath,[cfg.EEG_fileS '_DBS' num2str(DBS) '.sef']),data(:,Seg(i,1):Seg(i,2)),headerS);
            else
                header.dbs = false;
                NoDBS = NoDBS + 1;
                lab_write_sef(fullfile(cfg.DBS.EEG_filepath,[cfg.EEG_fileS '_NoDBS' num2str(NoDBS) '.sef']),data(:,Seg(i,1):Seg(i,2)),headerS);
            end
        end
    end
end

disp ('Finished')

end