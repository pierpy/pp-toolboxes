% Read Neuroscan cnt

function [data,header,cfg] = lab_read_cnt(filename)

[filenameL,filepath] = lab_filename(filename);

if ~check_header(filename,'Version 3.0')
    disp('    unkown file format')
    data = [];
    header = [];
    cfg = [];
end

orig = loadcnt(filename);
orig = rmfield(orig, {'data', 'ldnsamples'});

% do some reformatting/renaming of the header items
header.samplingrate = orig.header.rate;
header.numchannels = orig.header.nchannels;
header.numtimeframes = orig.header.nums;
for i=1:header.numchannels
    label{i} = deblank(orig.electloc(i).lab);
end
header.channels = char(label(1:header.numchannels));
clearvars label
header.version=0;
header.numauxchannels=0;
header.year=0;
header.month=0;
header.day=0;
header.hour=0;
header.minute=0;
header.second=0;
header.millisecond=0;
header.orig = orig;

tmp = loadcnt(filename,'sample1',0,'ldnsamples',header.numtimeframes);
data = tmp.data;

if isfield(orig,'event')
    events.POS = [];
    events.DUR = [];
    events.OFF = [];
    events.TYP = [];
    for i=1:numel(orig.event)
        events.TYP = [events.TYP cellstr('Marker')];
        events.POS = [events.POS int64(orig.event(i).offset + 1)];
        events.DUR = [events.DUR int64(0)];
        events.OFF = [events.OFF int64(0)];
        % the code above assumes that all events are stimulus markers
        % howevere, there are also interesting events possible, such as responses
        if orig.event(i).stimtype~=0
            events.TYP = [events.TYP cellstr(['Type' num2str(orig.event(i).stimtype)])];
            events.POS = [events.POS int64(orig.event(i).offset + 1)];
            events.DUR = [events.DUR int64(0)];
            events.OFF = [events.OFF int64(0)];
        end
        if orig.event(i).keypad_accept~=0
            events.TYP = [events.TYP cellstr(['Response' num2str(orig.event(i).keypad_accept)])];
            events.POS = [events.POS int64(orig.event(i).offset + 1)];
            events.DUR = [events.DUR int64(0)];
            events.OFF = [events.OFF int64(0)];
        end
    end
    header.events = events;
end
cfg.EEG_file = filenameL;
cfg.EEG_filepath= filepath;

return

function val = check_header(filename,head)
fid = fopen(filename, 'rb');
[str,nsiz] = fread(fid,length(head), 'uint8=>char');
fclose(fid);
if nsiz~=length(head)
    val = false;
else
    val = all(str(:)==head(:));
end
return