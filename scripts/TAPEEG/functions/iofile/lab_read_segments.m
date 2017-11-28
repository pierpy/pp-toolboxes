% written by F. Hatz 2016

function header = lab_read_segments(Filename,header)

if ~exist('header','var')
    header = [];
end
if ~exist('Filename','var')
    [Filename,Filepath] = uigetfile('','Select EEG-File','*.sef');
    Filename = fullfile(Filepath,Filename);
end
disp(['    read info from ' Filename])
[~,Filepath,~,FilenameS] = lab_filename(Filename);

fid=fopen(fullfile(Filepath,[FilenameS '_segment.vrb']),'r');
if fid > 0
    tline=fgetl(fid);
    if length(tline) < 11 | ~strcmp(tline(1:11),'Segment EEG')
        return
    end
    tline=fgetl(fid);
    events.POS = [];
    events.DUR = [];
    events.OFF = [];
    events.TYP = {};
    while ~isnumeric(tline)
        if length(tline) > 6 & strcmp(tline(1:7),'Segment')
             tmp = textscan(tline,'%s');
             tmp = tmp{1};
             if size(tmp,1) == 5
                 events.POS(1,end+1) = int64(str2num(tmp{3})); %#ok<ST2NM>
                 events.DUR(1,end+1) = int64(1);
                 events.OFF(1,end+1) = int64(0);
                 events.TYP{1,end+1} = ['START_' tmp{2}(1:end-1)];
                 events.POS(1,end+1) = int64(str2num(tmp{5})); %#ok<ST2NM>
                 events.DUR(1,end+1) = int64(1);
                 events.OFF(1,end+1) = int64(0);
                 events.TYP{1,end+1} = ['STOP_' tmp{2}(1:end-1)];
             end
        end
        tline=fgetl(fid);
    end
    if ~isempty(events.POS)
        events.OFF = int64(zeros(1,size(events.POS,2)));
        events.POS = int64(events.POS);
        events.DUR = int64(events.DUR);
    end
    header = lab_mix_markers(header,events);
end