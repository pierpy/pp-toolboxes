% function to stitch to EEG/MEG files
%
% [data,header] = lab_stitch(data1,header1,data2,header2,window)
%
% data1 & header1       first EEG/MEG (Output of lab_read_data)
% data2 & header2       first EEG/MEG (Output of lab_read_data)
% window                if single value: hanning window of .. seconds is
%                                        created for stitching
%                       if vector: values must be between 0 and 1
%
% written by F. Hatz, 2013

function [data,header] = lab_stitch(data1,header1,data2,header2,window)

if header1.numchannels ~= header2.numchannels
    data = data1;
    header = header1;
    disp('    Stitching not possible, different number of channels')
    return
end

if ~exist('window','var')
    disp('    Stitching with 1 second hanning window')
    window = 1-hann(header1.samplingrate);
elseif length(window) == 1 & isnumeric(window) & window > 0
    disp(['    Stitching with ' num2str(window) ' second hanning window'])
    window = 1-hann(ceil(window*header1.samplingrate));
elseif length(window) == 1 & isnumeric(window) & window <= 0
    % disable window
    disp('    Stitching without window')
    window = ones(1,4);
end
if size(window,2) == 1
    window = window';
end
if size(window,1) ~= 1
    data = data1;
    header = header1;
    disp('    Stitching not possible, wrong window information')
    return
end
if mod(size(window,2),2) > 0
    window = [window window(end)];
end

% Correct saplingrate of second dataset if necessary
if header1.samplingrate ~= header2.samplingrate
     [data2,header2]=lab_resample_data(data2,header2,header1.samplingrate);
end

% Stitch data
Lwindow = length(window);
data = [data1 data2];
Sstop = size(data1,2) - Lwindow/2;
Sstart = Lwindow/2 + 1;
tmp = data1(:,Sstop+1:end) .* repmat(window(Lwindow/2+1:end),[size(data,1) 1]) + ...
    data2(:,1:Sstart-1) .* repmat(window(1:Lwindow/2),[size(data,1) 1]);
data = [data1(:,1:Sstop) tmp data2(:,Sstart:end)];

% Create final header / mix markers
header = header1;
header.numtimeframes = size(data,2);
if isfield(header2,'events') & isfield(header2.events,'POS') & ~isempty(header2.events.POS)
    header2.events.POS = header2.events.POS + size(data1,2) - Lwindow/2;
    header = lab_mix_markers(header,header2);
end
if isfield(header,'badchans') & isfield(header2,'badchans')
    header.badchans = union(header.badchans,header2.badchans);
    if size(header.badchans,1) > 1
        header.badchans = header.badchans';
    end
    header.goodchans = setdiff(header.goodchans,header.badchans);
    header.goodchans = header.goodchans(:)';
end
if isfield(header,'activationsexcluded') & isfield(header2,'activationsexcluded')
    header.activationsexcluded = max(header.activationsexcluded,header2.activationsexcluded);
end

% Create marker for stitch region
headertmp.events.POS = int64(Sstop+1);
headertmp.events.DUR = int64(Lwindow);
headertmp.events.OFF = int64(0);
headertmp.events.TYP = cellstr('Stitch');
header = lab_mix_markers(header,headertmp);

return
