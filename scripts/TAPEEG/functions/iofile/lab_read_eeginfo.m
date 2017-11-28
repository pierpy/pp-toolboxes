% header = lab_read_eeginfo(Filename,header) - read TAPEEG eeginfo *.info
% written by F. Hatz 2012

function header = lab_read_eeginfo(Filename,header,novrb)

if ~exist('novrb','var')
    novrb = false;
end
if ~exist('header','var')
    header = [];
end

if novrb  == false
    disp(['    read info from ' Filename])
end
[~,Filepath,~,FilenameS] = lab_filename(Filename);

fid=fopen(fullfile(Filepath,[FilenameS '.info']),'r');
if fid > 0
    tline=fgetl(fid);
    while ~isnumeric(tline)
        switch tline
            case 'Bad channels'
                header.badchans = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Bad activations'
                header.badchans = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Good channels'
                header.goodchans = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Interpolated channels'
                header.interpolated = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Reference channel'
                tmp = fgetl(fid);
                if ~isnan(str2double(tmp))
                    header.ref_chan = str2num(tmp); %#ok<ST2NM>
                elseif ~isempty(tmp)
                    header.ref_chan = tmp;
                else
                    header.ref_chan = [];
                end
            case 'Samplingrate'
                if ~isfield(header,'samplingrate') | isempty(header.samplingrate)
                    header.samplingrate = str2num(fgetl(fid)); %#ok<ST2NM>
                end
            case 'Highpass filter'
                header.highpass = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Lowpass filter'
                header.lowpass = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'ECG channel'
                header.ecg_ch = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Number EEG channels'
                if ~isfield(header,'numchannels') | isempty(header.numchannels)
                    header.numchannels = str2num(fgetl(fid)); %#ok<ST2NM>
                end
            case 'Number AUX channels'
                header.numauxchannels = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Datatype'
                header.datatype = fgetl(fid);
            case 'Starttime in seconds'
                header.begin = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Stoptime in seconds'
                header.end = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Excluded activations'
                header.activationsexcluded = str2num(fgetl(fid)); %#ok<ST2NM>
            case 'Excluded channels:'
                header.badchans = str2num(fgetl(fid)); %#ok<ST2NM>
                if isfield(header,'numchannels')
                    header.goodchans = setdiff(1:header.numchannels,header.badchans);
                    header.goodchans = header.goodchans(:)';
                    header.ref_chan = header.numchannels;
                end
            case 'Quality'
                header.quality = str2num(fgetl(fid)); %#ok<ST2NM>
        end
        tline=fgetl(fid);
    end
    if isfield(header,'begin') & isempty(header.begin)
        header = rmfield(header,'begin');
    end
    if isfield(header,'end') & isempty(header.end)
        header = rmfield(header,'end');
    end
    if isfield(header,'activationsexcluded') & isempty(header.activationsexcluded)
        header = rmfield(header,'activationsexcluded');
    end
    fclose(fid);
    clearvars tline fid
    
    if isfield(header,'badchans') & max(header.badchans == 0);
        if ~isempty(find(header.badchans ~= 0,1))
            header.badchans = header.badchans(header.badchans ~= 0);
        else
            header.badchans = [];
        end
    end
    if isfield(header,'badchans') & isfield(header,'numdatachannels')
        header.goodchans = setdiff(1:header.numdatachannels,header.badchans);
        header.goodchans = header.goodchans(:)';
    end
    if isfield(header,'datatype') & header.datatype == -1
        header.datatype = 'eeg';
    end
    if isfield(header,'numauxchannels') & isfield(header,'numdatachannels') & header.numauxchannels + header.numdatachannels ~= header.numchannels
        header.numdatachannels = header.numchannels - header.numauxchannels;
    end
end