% Save eeg/meg data to different file formats
%
% [filename,settings] = lab_save_data(data,header,format,filename,settings,cfg)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% format   = extension of fileformat to store data (optional)
% filename = path and filename to store data (optional)
% settings = structure with config (optional)
% cfg      = structure with config (optional)
%
% written by F. Hatz 2012

function [filename,settings] = lab_save_data(data,header,format,filename,settings,cfg)

% check if only filename and no format was entered
if exist('format','var')
    tmp = strfind(format,filesep);
    if ~isempty(tmp)
        filename = format;
        clearvars format
    else
        tmp = strfind(format,'.');
        if ~isempty(tmp) & tmp(1)>1
            filename = format;
            clearvars format
        elseif ~isempty(tmp) & tmp == 1 & length(format) > 1
            format = format(2:end);
        end
    end
    clearvars tmp
end

if ~exist('filename','var')
    if exist('format','var')
        [filename,filepath]=uiputfile(['*.' format],'Select file and path');
    else
        [filename,filepath]=uiputfile('*.*','Select file and path');
    end
    filename = fullfile(filepath,filename);
    clearvars filepath;
end
if ~exist('format','var')
    [~,~,format] = lab_filename(filename);
end

% correct for missing header
if ~exist('header','var')
    if isnumeric(data)
        header.numchannels = size(data,1);
        header.numdatachannels = size(data,1);
        header.numauxchannels = 0;
        header.ref_chan = 'unkown';
        header.numtimeframes = size(data,2);
        header.samplingrate = [];
        header.channels = num2str((1:header.numchannels)');
        header.datatype = 'unkown';
    else
        header = [];
    end
end
if ~exist('settings','var')
    settings = [];
end

% correct format
format = lower(format);
tmp = strfind(format,'.');
if ~isempty(tmp) & tmp(end) ~= length(format)
    format = format(tmp(end)+1:end);
end
clearvars tmp
if strcmp(format,'fiff')
    if isfield(header,'orig') & isfield(header.orig,'info')
        testfiff = 1;
    else
        testfiff = 0;
    end
    if testfiff == 0
        disp('  No FIFF-header-info found, writing sef-file')
        format = 'sef';
    end
    clearvars testfiff
end

% correct filename
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.' format]);
clearvars tmp filepath

if exist(['lab_write_' format]) == 2 %#ok<EXIST>
    % check for scale factor if necessary
    if strcmp(format,'txt')
        if ~isfield(settings,'scaletxt')
            if exist('cfg','var') & isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
                settings.scaletxt = 'auto';
            else
                disp ('   Ask for scale factor')
                prompt={'Scale data for txt-export (auto / number / 1 = off)'};
                name='Scale factor txt export';
                numlines(1,1) = 1;
                numlines(1,2) = 50;
                defaultanswer={'1'};
                answer=inputdlg(prompt,name,numlines,defaultanswer);
                pause(0.1)
                if size(answer,1) > 0
                    settings.scaletxt = answer{1,1};
                else
                    settings.scaletxt = '1';
                end
                clearvars defaultanswer answer name numlines prompt
            end
        end
        if ischar(settings.scaletxt) & strcmp(settings.scaletxt,'auto')
            settings.scaletxt = round(1000/max(max(abs(data))));
        elseif ischar(settings.scaletxt) & ~isnan(str2double(settings.scaletxt)) & str2num(settings.scaletxt) > 0 %#ok<ST2NM>
            settings.scaletxt = str2num(settings.scaletxt); %#ok<ST2NM>
        elseif ~isnumeric(settings.scaletxt) | isempty(settings.scaletxt)
            settings.scaletxt = 1;
        elseif isnumeric(settings.scaletxt) & settings.scaletxt <= 0
            settings.scaletxt = 1;
        end
        if settings.scaletxt > 10
            data = round(data * settings.scaletxt); %#ok<NASGU>
            format = 'txt_int';
        else
            data = data * settings.scaletxt; %#ok<NASGU>
        end
    end
    
    % script for writing files is selected by name (all scripts in folder 'iofile')
    if ~isempty(header)
        eval(['lab_write_' format '(filename,data,header);']);
    else
        eval(['lab_write_' format '(filename,data);']);
    end
else
    disp('Error writing: data-format not supported')
end


