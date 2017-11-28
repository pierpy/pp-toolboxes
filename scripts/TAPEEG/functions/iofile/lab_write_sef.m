% Write eeg/meg to Cartool .sef
%
% lab_write_sef(filename,data,header)
%
% written by F. Hatz 2012

function lab_write_sef(filename,data,header)

if ~exist('header','var')
    header = lab_create_header(data);
end

[~,filepath,~,filenameS] = lab_filename(filename);
filename = [filenameS '.sef'];
header.version = 'SE01';

if ~isfield(header,'year')
    header.year = 0;
    header.month = 0;
    header.day = 0;
    header.hour = 0;
    header.minute = 0;
    header.second = 0;
    header.millisecond = 0;
end

fid=fopen(fullfile(filepath,filename),'w');

% Write header
fwrite(fid,header.version,'int8');
fwrite(fid,header.numchannels,'int32');
fwrite(fid,header.numauxchannels,'int32');
fwrite(fid,header.numtimeframes,'int32');
fwrite(fid,header.samplingrate,'float32');
fwrite(fid,header.year,'int16');
fwrite(fid,header.month,'int16');
fwrite(fid,header.day,'int16');
fwrite(fid,header.hour,'int16');
fwrite(fid,header.minute,'int16');
fwrite(fid,header.second,'int16');
fwrite(fid,header.millisecond,'int16');

% Write channel names
channels = zeros(header.numchannels,8);
for i=1:header.numchannels
    tmp=uint8(header.channels(i,:));
    if size(tmp,2) <= 8
        channels(i,1:size(tmp,2)) = tmp;
    else
        channels(i,:) = tmp(1,1:8);
    end
end
fwrite(fid,channels','int8');

% Write data
fwrite(fid,data,'float32');
fclose(fid);

% Write Marker-file
if isfield(header,'events')
    Marker_file=fullfile(filepath,[filename '.mrk']);
    lab_write_mrk(Marker_file,header);
end

if isfield(header,'ref_chan') 
    % Write EEGinfo-file (*.txt)
    if length(filename) <= 10 | ~strcmp(filename(end-10:end),'ICAtopo.sef')
        lab_write_eeginfo(fullfile(filepath,filename),header)
    end
end

% Write loc file
ELS_file = fullfile(filepath,[filenameS '.els']);
ELS_file = lab_write_locs(ELS_file,header,'bad');
if ~isempty(ELS_file)
    % Write *.LM-file
    fidout=fopen(fullfile(filepath,[filenameS '.lm']),'w');
    fprintf(fidout,[filename native2unicode([13 10])]);
    fprintf(fidout,[ELS_file native2unicode([13 10])]);
    fclose(fidout);
end

% Write weights
if isfield(header,'W')
    dlmwrite(fullfile(filepath,[filenameS '.w']),header.W,'delimiter','\t','precision', 6);
end

% Write individual freqbands
if isfield(header,'IFREQ') & ~isempty(header.IFREQ)
    IFREQ = header.IFREQ; %#ok<NASGU>
    save(fullfile(filepath,[filenameS '.ifreq']),'IFREQ');
end