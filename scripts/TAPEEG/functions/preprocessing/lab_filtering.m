function [data,header,cfg,skipprocessing] = lab_filtering(data,header,cfg,novrb)

skipprocessing = 0;
    
if ~exist('novrb','var')
    novrb = 0;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('header','var')
    tmp = inputdlg({'Samplingrate'},'Samplingrate',[1 25]);
    header.samplingrate = str2num(tmp{1}); %#ok<ST2NM>
    header.numchannels = size(data,1);
    header.numdatachannels = size(data,1);
    clearvars tmp
    novrb = 1;
elseif isnumeric(header)
    tmp.samplingrate = header;
    tmp.numchannels = size(data,1);
    tmp.numdatachannels = size(data,1);
    header = tmp;
    clearvars tmp
    novrb = 1;
end

if novrb == 0
    if ~isfield(cfg,'EEG_file') & isfield(header,'EEG_file')
        cfg.EEG_file = header.EEG_file;
        cfg.EEG_filepath = header.EEG_filepath;
    end
    if isfield(cfg,'EEG_file')
        [~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
    else
        novrb = 1;
    end
end

if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end

if ~isfield(cfg,'FILT') | isempty(cfg.FILT)
    [cfg,skipprocessing] = lab_set_filtering(cfg,header);
    if skipprocessing ~= 0
        return
    end
end

if ~isfield(cfg.FILT,'Filter') | isempty(cfg.FILT.Filter)
    disp('Skip Filtering - invalid settings')
    return
end

% Write verbose file
if novrb == 0
    Verbose_file=[cfg.EEG_fileS '_filt.vrb'];
    fid=fopen(fullfile(cfg.EEG_filepath,Verbose_file),'w');
    fprintf(fid,'Filter EEG\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n');
    fprintf(fid,['EEG input file: ' cfg.EEG_file]);
    fprintf(fid,'\n\n');
    fprintf(fid,['Filtermode: ' cfg.FILT.filtermode]);
    fprintf(fid,'\n\n');
end

settings = cfg.FILT;
for i = 1:length(cfg.FILT.Filter)
    if ~isempty(cfg.FILT.Filter(i).firstchan)
        FirstChan = cfg.FILT.Filter(i).firstchan;
        if ~isempty(cfg.FILT.Filter(i).lastchan)
            LastChan = cfg.FILT.Filter(i).lastchan;
        else
            LastChan = FirstChan;
        end
        if LastChan > size(data,1)
            disp(['   Last channel not in data-range, set last channel to ' num2str(size(data,1))])
            LastChan = size(data,1);
        end
        
        if FirstChan == round(FirstChan) & FirstChan >= 1 & FirstChan <= size(data,1) & ...
                LastChan == round(LastChan) & LastChan >= 1
            datatmp = data(FirstChan:LastChan,:);
            headertmp = header;
            headertmp.numchannels = size(datatmp,1);
            headertmp.numdatachannels = headertmp.numchannels;
            
            settings.notch = cfg.FILT.Filter(i).notch;
            settings.highpass = cfg.FILT.Filter(i).highpass;
            settings.lowpass = cfg.FILT.Filter(i).lowpass;
            
            [datatmp,headertmp] = lab_filter(datatmp,headertmp,settings,1);
            
            data(FirstChan:LastChan,:) = datatmp;
            if isfield(headertmp,'notch')
                if ~isfield(header,'notch')
                    header.notch = [];
                end
                header.notch = [header.notch headertmp.notch];
            end
            if isfield(headertmp,'highpass')
                header.highpass = headertmp.highpass;
                header.lowpass = headertmp.lowpass;
            end
            if isfield(headertmp,'filter')
                header.filter = headertmp.filter;
            end
            
            if novrb == 0
                fprintf(fid,['From ' num2str(FirstChan) ' Channel to ' num2str(LastChan) ' Channel:']);
                fprintf(fid,'\n');
                fprintf(fid,'Notch Filter: ');
                for j = 1:size(settings.notch,2)
                    fprintf(fid,[num2str(settings.notch(1,j)) 'Hz ']);
                end
                fprintf(fid,'\n');
                fprintf(fid,['Highpass Filter: ' num2str(settings.highpass) 'Hz']);
                fprintf(fid,'\n');
                fprintf(fid,['Lowpass Filter: ' num2str(settings.lowpass) 'Hz']);
                fprintf(fid,'\n\n');
            end
        else
            disp('Skip Filtering - invalid range of channels')
            if novrb == 0
                fprintf(fid,['Filtering from ' num2str(FirstChan) ' Channel to ' num2str(LastChan) ' Channel not possible']);
            end
        end
    end
end

if novrb == 0
    if isfield(header,'filter') & isfield(header.filter,'order')
        fprintf(fid,'\n\n');
        fprintf(fid,['Filterorder: ' num2str(header.filter.order)]);
    elseif isfield(header,'filter') & isfield(header.filter,'smoothingfactor')
        fprintf(fid,'\n\n');
        fprintf(fid,['Smoothing factor: ' num2str(header.filter.smoothingfactor)]);
    end
    fclose(fid);
end

end