% Read Cartool mrk-files
%
% header = lab_read_mrk(filename,header)
%
% header   = output of lab_read_data (if variable header is defined, marker
%            info stored in header will be combined with new marker info 
%            of loaded file)
%
% To shift timings in the mrk-file to match an altered eeg/meg-file, you
% can add in the first line the number of timeframes to shift (positive and
% negative values possible)
%                TL02 \t shift (\t = tabulator)
%
% written by F. Hatz 2012

function header = lab_read_mrk(filename,header)

if ~exist('header','var') | ~isfield(header,'numtimeframes')
    header.numtimeframes = 1000000000000;
    eventsonly = 1;
else
    eventsonly = 0;
end
events.POS = [];
events.DUR = [];
events.OFF = [];
events.TYP = {};

fid=fopen(filename,'r');
if fid > 0
    tline=fgets(fid);
    if strcmp(tline(1:4),'TL02')
        if ~isnumeric(tline)
            tmp = regexp(tline,'\t', 'once');
            if ~isempty(tmp)
                shift = str2num(tline(tmp(1)+1:end)); %#ok<ST2NM>
                disp(['    correct marker positions for detected skip ' num2str(shift)])
            else
                shift = 0;
            end
            tline=fgets(fid);
        end
        while ~isnumeric(tline)
            if regexp(tline,'\t','once')
                tmp=regexp(tline,'\t');
                if size(tmp,2) == 3
                    position = str2num(tline(tmp(1)+1:tmp(2)-1)) - shift + 1; %#ok<ST2NM>
                    duration = str2num(tline(tmp(2)+1:tmp(3)-1)) - shift + 1; %#ok<ST2NM>
                    if position <= header.numtimeframes & position > 0;
                        events.POS(1,end+1) = int64(position);
                        events.DUR(1,end+1) = int64(duration + 1) - events.POS(1,end);
                        tmp = regexp(tline,'"');
                        events.TYP(1,end+1) = cellstr(tline(tmp(1)+1:tmp(2)-1));
                    end
                    clearvars position duration
                else
                    position = str2num(tline(1:tmp(1)-1)) - shift + 1; %#ok<ST2NM>
                    duration = str2num(tline(tmp(1)+1:tmp(2)-1)) - shift + 1; %#ok<ST2NM>
                    if position <= header.numtimeframes & position > 0;
                        events.POS(1,end+1) = int64(position);
                        events.DUR(1,end+1) = int64(duration) - events.POS(1,end);
                        tmp = regexp(tline,'"');
                        events.TYP(1,end+1) = cellstr(tline(tmp(1)+1:tmp(2)-1));
                    end
                    clearvars position duration
                end
                clearvars tmp
            end
            tline=fgets(fid);
        end
        fclose(fid);
        clearvars tline tmp i
        if ~isempty(events.POS)
            events.OFF = int64(zeros(1,size(events.POS,2)));
            events.POS = int64(events.POS);
            events.DUR = int64(events.DUR);
        end
    else
        disp('    error: mrk-fileformat not supported')
    end
end
clearvars fid;

if eventsonly == 1
    header = events;
elseif ~isempty(events)
    if isfield(header,'events') & ~isempty(header.events)
        header = lab_mix_markers(header,events);
    else
        header.events = events;
    end
end