function [header] = lab_reduce_header(header,includechannels)

if isfield(header,'goodchans')
    % redefine goodchans & badchans
    allchans=[1:header.numchannels;zeros(1,header.numchannels)];
    allchans(2,header.goodchans)=1;
    if isfield(header,'interpolated')
        allchans(2,header.interpolated)=2;
    end
    if isfield(header,'badchans') & ~isempty(header.badchans) & min(header.badchans) ~= 0
        allchans(2,header.badchans)=3;
    end
    allchans = allchans(:,includechannels);
    allchans(1,:) = 1:size(allchans,2);
    header.goodchans=find(allchans(2,:)==1);
    if isfield(header,'interpolated')
        header.interpolated = find(allchans(2,:)==2);
    end
    header.badchans=find(allchans(2,:)==3);
end
if isfield(header,'bad') & isfield(header.bad,'epochs')
    header.bad.epochs = header.bad.epochs(includechannels,:);
end

% redefine ref_chan & ecg_ch
if isnumeric(header.ref_chan)
    if header.ref_chan > 0
        for i = 1:size(header.ref_chan,2)
            if header.ref_chan(1,i) > 0
                tmp = find(includechannels == header.ref_chan(1,i));
                if ~isempty(tmp)
                    header.ref_chan(1,i) = tmp;
                end
                clearvars tmp
            end
        end
    end
end
if isfield(header,'ecg_ch') & header.ecg_ch > 0
    header.ecg_ch = find(includechannels == header.ecg_ch);
    if isempty(header.ecg_ch)
        header.ecg_ch = 0;
    end
end
if isfield(header,'eog_ch') & header.eog_ch > 0
    header.eog_ch = find(includechannels == header.eog_ch);
    if isempty(header.eog_ch)
        header.eog_ch = 0;
    end
end

% redefine channel labels
header.channels = header.channels(includechannels,:);

% Correct header
if isfield(header,'includechans')
    header.includechans = header.includechans(1,includechannels);
else
    header.includechans = includechannels;
end
if max(includechannels) <= header.numdatachannels | header.numdatachannels == header.numchannels
    header.numchannels = length(includechannels);
    header.numdatachannels = header.numchannels;
    header.numauxchannels = 0;
    header.includedatachans = header.includechans;
else
    header.numchannels = length(includechannels);
    header.numauxchannels = length(find(includechannels > header.numdatachannels));
    header.numdatachannels = header.numchannels - header.numauxchannels;
    header.includedatachans = header.includechans(1,1:end-header.numauxchannels);
end
if isfield(header,'ch_names')
    header.ch_names = header.ch_names(includechannels,1);
end

% redefine locs
if isfield(header,'locs')
    includedatachannels = includechannels(1,1:end-header.numauxchannels);
    header.locs.x = header.locs.x(1,includedatachannels);
    header.locs.y = header.locs.y(1,includedatachannels);
    header.locs.z = header.locs.z(1,includedatachannels);
    header.locs.labels = header.locs.labels(1,includedatachannels);
    if isfield(header.locs,'radius')
        header.locs.radius = header.locs.radius(1,includedatachannels);
        header.locs.theta = header.locs.theta(1,includedatachannels);
    end
end

return