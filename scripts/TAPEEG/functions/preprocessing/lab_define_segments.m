% Split EEG/MEG in segments
%
% [Segments,cfg] = lab_define_segments(data,header,cfg)
%
% data     = eeg/meg data (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
%
% Written by F. Hatz 2012 Neurology Basel

function [Segments,cfg,skipprocessing] = lab_define_segments(data,header,cfg)

global DATA_SEG HEADER_SEG SEGAUTO
    
disp ('Define segments for analysis')
skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end
if ~isfield(cfg,'settings_path')
    cfg.settings_path = pwd;
end
[~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);

if ~exist('cfg','var') | ~isfield(cfg,'SEG') | isempty(cfg.SEG)
    [cfg,skipprocessing] = lab_set_define_segments(cfg,header);
    pause(0.2);
    if skipprocessing == 1
        Segments = [];
        return
    end
end

% convert old settings
if isnumeric(cfg.SEG.select)
    tmp = {'Select complete file';'Select by Markers';'Select by Start-Stop-Marker';'Select by EDF-file';'Automated'};
    cfg.SEG.select = tmp{cfg.SEG.select};
    clearvars tmp;
end

if strcmp(cfg.SEG.select,'Select complete file');
    header.segments = [1 header.numtimeframes];
elseif strcmp(cfg.SEG.select,'Select by Markers') | strcmp(cfg.SEG.select,'Select by Start-Stop-Marker') | ...
        strcmp(cfg.SEG.select,'Select by EDF-file')
    if strcmp(cfg.SEG.select,'Select by EDF-file') & ~isempty(data)
        if ~isfield(header,'events') | ~isfield(header.events,'edfimport')
            disp('   Missing *.edf')
            if ~exist(fullfile(cfg.EEG_filepath,[cfg.EEG_fileS '~.edf']),'file')
                disp('   Write *~.edf')
                cfg.SEG.EDF.interpolatebad = true;
                cfg.SEG.EDF.filepath = '';
                cfg.SEG.EDF.mountpath = '';
                cfg.SEG.EDF_file=[cfg.EEG_fileS '~.edf'];
                cfg.SEG.EDF_filepath = cfg.EEG_filepath;
                [~,~,cfg.SEG] = lab_export2edf(data,header,cfg.SEG,'nofolder');
                disp(['   Please edit ' cfg.SEG.EDF_file ', rename to ' cfg.EEG_fileS '.edf and restart'])
                clearvars cfgtmp
            end
            Segments = [];
            return
        end
    end
    novalid = 0;
    if strcmp(cfg.SEG.select,'Select by Markers')
        if ~isfield(cfg.SEG,'markerstart') | isempty(cfg.SEG.markerstart)
            novalid = 1;
        end
    else
        if ~isfield(cfg.SEG,'markerstart') | isempty(cfg.SEG.markerstart) | ...
                ~isfield(cfg.SEG,'markerstop') | isempty(cfg.SEG.markerstop)
            novalid = 1;
        end
    end
    if novalid == 1
        if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 0
            disp('   Abort: No valid markers for segments selected')
            Segments = [];
            return
        else
            if ~exist(fullfile(cfg.settings_path,'segments.xls'),'file')
                tmp = unique(header.events.TYP);
                markerselection = [cellstr('*start*') tmp(:)' cellstr('*end*')];
                lab_write_xls(fullfile(cfg.settings_path,'segments~.xls'),markerselection);
                clearvars markerselection tmp
            end
            while ~exist(fullfile(cfg.settings_path,'segments.xls'),'file')
                disp('   Please select start&end marker in -segments~.xls- and rename to -segments.xls-')
                pause(10);
            end
            if ispc
                [~,xlsinput] = xlsread(fullfile(cfg.settings_path,'segments.xls'));
            else
                [~,xlsinput] = xlsread(fullfile(cfg.settings_path,'segments.xls'),1,'','basic');
            end
            cfg.SEG.markerstart = xlsinput{1,1};
            cfg.SEG.markerstop = xlsinput{1,2};
        end
    end
    if isfield(header,'events') & ~isempty(header.events)
        if strcmp(cfg.SEG.select,'Select by Markers')
            segments = [];
            if isfield(cfg.SEG,'markerstart') & ~isempty(cfg.SEG.markerstart)
                markerstart = header.events.POS(1,ismember(header.events.TYP,cfg.SEG.markerstart)==1);
                tmp = header.events.DUR(1,ismember(header.events.TYP,cfg.SEG.markerstart)==1);
                if min(tmp) > 1
                    markerend = markerstart + tmp;
                    segments = cat(1,segments,[markerstart' markerend']);
                end
                clearvars markerstart markerend tmp
            end
            if isfield(cfg.SEG,'markerstop') & ~isempty(cfg.SEG.markerstop)
                markerstart = header.events.POS(1,ismember(header.events.TYP,cfg.SEG.markerstop)==1);
                tmp = header.events.DUR(1,ismember(header.events.TYP,cfg.SEG.markerstop)==1);
                if min(tmp) > 1
                    markerend = markerstart + tmp;
                    segments = cat(1,segments,[markerstart' markerend']);
                end
                clearvars markerstart markerend tmp
            end
            if ~isempty(segments)
                header.segments = segments;
            else
                header.segments = [1 header.numtimeframes];
                disp('   No valid markers, process complete file')
            end
        else
            if isfield(cfg.SEG,'markerstart') & ~isempty(cfg.SEG.markerstart) & ...
                    isfield(cfg.SEG,'markerstop') & ~isempty(cfg.SEG.markerstop)
                if strcmp(cfg.SEG.markerstop,'*end*') & strcmp(cfg.SEG.markerstart,'*start*')
                    markerstart = 1;
                    markerend = header.numtimeframes;
                elseif strcmp(cfg.SEG.markerstart,'*start*')
                    markerend = header.events.POS(1,ismember(header.events.TYP,cfg.SEG.markerstop)==1);
                    markerstart = ones(1,size(markerend,2));
                elseif strcmp(cfg.SEG.markerstop,'*end*') & ~strcmp(cfg.SEG.markerstart,'*start*')
                    markerstart = header.events.POS(1,ismember(header.events.TYP,cfg.SEG.markerstart)==1);
                    markerend = ones(1,size(markerstart,2))*header.numtimeframes;
                else
                    markerstart = header.events.POS(1,ismember(header.events.TYP,cfg.SEG.markerstart)==1);
                    markerend = header.events.POS(1,ismember(header.events.TYP,cfg.SEG.markerstop)==1);
                end
                if isempty(markerstart)
                    markerstart = 1;
                    markerend = header.numtimeframes;
                end
                if size(markerstart,2) == size(markerend,2)
                    header.segments = [markerstart' markerend'];
                else
                    header.segments = [1 header.numtimeframes];
                    disp('   No equal number of start- and end-markers, process complete file')
                end
                clearvars markerstart markerend
            end
        end
    else
        header.segments = [1 header.numtimeframes];
        disp('   No markers in file, process complete file')
    end
    clearvars headertmp
elseif strcmp(cfg.SEG.select,'Automated')
    if isfield(cfg.SEG.AUTO,'mergefolder') & cfg.SEG.AUTO.mergefolder == true & ...
            isfield(cfg,'lastfilefolder') & isfield(cfg,'firstfilefolder')        
        if cfg.firstfilefolder == true
            DATA_SEG = [];
            HEADER_SEG = [];
        end
        if cfg.lastfilefolder == false
            if isempty(DATA_SEG)
                DATA_SEG = data;
                HEADER_SEG = header;
            else
                [DATA_SEG,HEADER_SEG] = lab_stitch(DATA_SEG,HEADER_SEG,data,header,1);
            end
            Segments = [];
            return
        else
            if ~isempty(DATA_SEG)
                [data,header] = lab_stitch(DATA_SEG,HEADER_SEG,data,header,1);
            end
            DATA_SEG = [];
            HEADER_SEG = [];
        end
    end
    header.segments = [];
    if ~isempty(SEGAUTO) & isfield(cfg,'EEG_filepath') & isfield(cfg,'EEG_fileS')
        tmp = find(strcmp(SEGAUTO(1,:),fullfile(cfg.EEG_filepath,cfg.EEG_fileS)));
        tmp2 = find(cell2mat(SEGAUTO(2,tmp)) == header.numtimeframes);
        if ~isempty(tmp) & ~isempty(tmp2)
            disp('   Take auto-segments from previous calculation');
            header.segments = SEGAUTO{3,tmp(tmp2(1))};
        end
        clearvars tmp tmp2
    end
    if isempty(header.segments)
        header.segments = lab_segments_auto(data,header,cfg);
        if isfield(cfg,'EEG_filepath') & isfield(cfg,'EEG_fileS')
            if isempty(SEGAUTO)
                SEGAUTO{1,1} = fullfile(cfg.EEG_filepath,cfg.EEG_fileS);
                SEGAUTO{2,1} = header.numtimeframes;
                SEGAUTO{3,1} = header.segments;
            else
                SEGAUTO{1,end+1} = fullfile(cfg.EEG_filepath,cfg.EEG_fileS);
                SEGAUTO{2,end} = header.numtimeframes;
                SEGAUTO{3,end} = header.segments;
            end
        end
    end
else
    header.segments = [];
end

if ~isempty(header.segments)
    header.segments = int64(header.segments);
    for i = 1:size(header.segments,1)
        disp(['   Create segment ' num2str(header.segments(i,1)) ' - ' num2str(header.segments(i,2))]);
        headertmp = header;
        if ~isempty(data) & header.segments(i,2) <= size(data,2)
            Segments{1,i} = data(:,(header.segments(i,1):header.segments(i,2))); %#ok<AGROW>
            headertmp.numtimeframes = size(Segments{1,i},2);
        else
            headertmp.numtimeframes = header.segments(i,2) - header.segments(i,1) + 1;
        end
        
        if isfield(header,'events') & isfield(header.events,'POS') & ...
                ~isempty(header.events.POS) & ~isempty(header.events.DUR)
            headertmp.events.POS =  header.events.POS - header.segments(i,1) + 1;
            tmp = find(headertmp.events.POS > 0);
            tmp2 = find(headertmp.events.POS <= headertmp.numtimeframes);
            tmp = intersect(tmp,tmp2);
            if ~isempty(tmp)
                headertmp.events.POS =  headertmp.events.POS(1,tmp);
                headertmp.events.TYP =  headertmp.events.TYP(1,tmp);
                headertmp.events.DUR =  headertmp.events.DUR(1,tmp);
                if ~isempty(headertmp.events.OFF)
                    headertmp.events.OFF =  headertmp.events.OFF(1,tmp);
                else
                    headertmp.events.OFF =  [];
                end
            else
                headertmp.events.POS =  [];
                headertmp.events.TYP =  [];
                headertmp.events.DUR =  [];
                headertmp.events.OFF =  [];
            end
        end
        Segments{2,i} = headertmp; %#ok<AGROW>
        clearvars headertmp tmp tmp2
    end
    
    if ~strcmp(cfg.SEG.select,'Select complete file')
        Verbose_file = fullfile(cfg.EEG_filepath,[cfg.EEG_fileS '_segment.vrb']);
        fid=fopen(Verbose_file,'w');
        fprintf(fid,'Segment EEG\n');
        fprintf(fid,datestr(now,0));
        fprintf(fid,'\n');
        if isfield(header,'EEG_file')
            fprintf(fid,['EEG input file: ' header.EEG_file]);
            fprintf(fid,'\n\n');
        end
        if isfield(cfg,'SEG') & isfield(cfg.SEG,'markerstart')
            fprintf(fid,['Marker start: ' cfg.SEG.markerstart]);
            fprintf(fid,'\n');
        end
        if isfield(cfg,'SEG') & isfield(cfg.SEG,'markerstop')
            fprintf(fid,['Marker end: ' cfg.SEG.markerstop]);
            fprintf(fid,'\n\n');
        end
        for i = 1:size(header.segments,1)
            fprintf(fid,['Segment ' num2str(i) ': ' num2str(header.segments(i,1)) ' - ' num2str(header.segments(i,2))]);
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
        fclose(fid);
    end
else
    Segments = [];
end

return