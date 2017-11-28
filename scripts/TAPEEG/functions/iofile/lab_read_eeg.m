% Reads .eeg Neuroscan & .eeg Brainvision

function [data,header,cfg] = lab_read_eeg(filename,cfg)

[filenameL,filepath,~,filenameS] = lab_filename(filename);

if check_header(filename,'Version 3.0')
    type = 'ns_eeg';
    % manufacturer = 'Neuroscan';
    % content = 'epoched EEG';
elseif exist(fullfile(filepath, [filenameS '.vhdr']), 'file')
    type = 'brainvision_eeg';
    % manufacturer = 'BrainProducts';
    % content = 'continuous EEG data';
else
    type = 'unkown';
end

if strcmp(type,'ns_eeg')
    orig = read_ns_hdr(filename);
    % do some reformatting/renaming of the header items
    header.samplingrate = orig.rate;
    header.numtimeframes = orig.npnt*orig.nsweeps;
    header.numchannels = orig.nchan;
    header.channels = char(orig.label(1:header.numchannels));
    header.version=0;
    header.numauxchannels=0;
    header.year=0;
    header.month=0;
    header.day=0;
    header.hour=0;
    header.minute=0;
    header.second=0;
    header.millisecond=0;
    % remember the original header details
    header.orig = orig;
    
    % read data
    tmp = read_ns_eeg(filename,orig.nsweeps);
    siz = [orig.nsweeps orig.nchan orig.npnt];
    data = reshape(tmp.data, siz); % ensure 3-D array
    data = permute(data,[2 3 1]);
    data = reshape(data,orig.nchan,orig.npnt*orig.nsweeps);
    
    % create events
    header.events.POS=[];
    header.events.DUR=[];
    header.events.OFF=[];
    header.events.TYP=[];
    for i=1:hdr.nTrials
      % the *.eeg file has a fixed marker value for each trial
      % furthermore each trial has additional fields like accept, correct, response and rt
      tmp = read_ns_eeg(filename,i);
      % create an event with the marker value
      header.events.TYP = [header.events.OFF cellstr(['trial' num2str(i)])];
      header.events.POS = [header.events.POS int64((i-1)*orig.npnt + 1)];
      header.events.OFF = [header.events.OFF int64(0)];
      header.events.DUR = [header.events.DUR int64(orig.npnt)];
      % create an event with the boolean accept/reject code
      header.events.TYP = [header.events.OFF cellstr(['accept' num2str(tmp.sweep.accept)])];
      header.events.POS = [header.events.POS int64((i-1)*orig.npnt + 1)];
      header.events.OFF = [header.events.OFF int64(0)];
      header.events.DUR = [header.events.DUR int64(orig.npnt)];
      % create an event with the boolean correct/false code
      header.events.TYP = [header.events.OFF cellstr(['correct' num2str(tmp.sweep.correct)])];
      header.events.POS = [header.events.POS int64((i-1)*orig.npnt + 1)];
      header.events.OFF = [header.events.OFF int64(0)];
      header.events.DUR = [header.events.DUR int64(orig.npnt)];
      % create an event with the response time
      header.events.TYP = [header.events.OFF cellstr(['response' num2str(i)])];
      header.events.POS = [header.events.POS int64((i-1)*orig.npnt + 1 + tmp.sweep.rt*orig.rate)];
      header.events.OFF = [header.events.OFF int64(0)];
      header.events.DUR = [header.events.DUR int64(orig.npnt)];
    end
    cfg.EEG_file = filenameL;
    cfg.EEG_filepath= filepath;
elseif strcmp(type,'brainvision_eeg')
    warning off %#ok<WNOFF>
    orig = read_brainvision_vhdr(fullfile(filepath, [filenameS '.vhdr']));
    warning on %#ok<WNON>
    header.samplingrate = orig.Fs;
    header.numchannels = orig.NumberOfChannels;
    header.channels = char(orig.label(1:header.numchannels));
    header.numtimeframes = orig.nSamples*orig.nTrials;
    header.orig        = orig;
    data = read_brainvision_eeg(filename,orig,1,header.numtimeframes,1:header.numchannels);
    if exist(fullfile(filepath,[filenameS '.vmrk']),'file')
        header.events = lab_read_vmrk(fullfile(filepath,[filenameS '.vmrk']));
    end
    cfg.EEG_file = filenameL;
    cfg.EEG_filepath= filepath;
else
    disp('    unkown file format')
    data = [];
    header = [];
    cfg = [];
end


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