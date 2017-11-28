% Reduce number of channels
%
% [data,header] = lab_reduce_channels(data,header,includechannels)
%
% data            = matrix (chans x timeframes)
% header          = output of lab_read_data
% includechannels = indices of channels to include
%
% written by F. Hatz 2012

function [data,header] = lab_reduce_channels(data,header,includechannels,skiprefwarn)

if ~exist('skiprefwarn','var')
    skiprefwarn = false; %#ok<NASGU>
end

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('includechannels','var')
    if isfield(header,'channels')
        strlist = cellstr(header.channels);
    else
        strlist = cellstr(num2str((1:size(data,1))'))';
    end
    includechannels = listdlg('PromptString','Select included channels:','SelectionMode','multiple', ...
    'ListString',strlist,'InitialValue',1:length(strlist));
end
if isempty(includechannels)
    return
end

if isfield(header,'goodchans')
    % redefine goodchans & badchans
    allchans=[1:header.numchannels;zeros(1,header.numchannels)];
    allchans(2,header.goodchans)=1;
    if isfield(header,'interpolated')
        allchans(2,header.interpolated)=2;
    end
    if isfield(header,'badchans') & ~isempty(header.badchans)
        allchans(2,header.badchans)=3;
    end
    allchans = allchans(:,includechannels);
    allchans(1,:) = 1:size(allchans,2);
    header.goodchans = find(allchans(2,:)==1);
    if isfield(header,'interpolated')
        header.interpolated = find(allchans(2,:)==2);
    end
    header.badchans = find(allchans(2,:)==3);
end
if isfield(header,'bad') & isfield(header.bad,'epochs')
    tmp = includechannels(includechannels<=size(header.bad.epochs,1));
    header.bad.epochs = header.bad.epochs(tmp,:);
    clearvars tmp
end

% redefine ref_chan & ecg_ch & eog_ch
if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
    if header.ref_chan > 0
        % if length(union(includechannels,header.ref_chan)) > length(includechannels) & ...
        %         skiprefwarn == false & max(header.ref_chan) <= header.numchannels
        %     disp('    warning: reference channel can not be omitted and is added')
        %     includechannels = union(includechannels,header.ref_chan);
        % end
        ref_chan = [];
        for i = 1:size(header.ref_chan,2)
            if header.ref_chan(1,i) > 0
                tmp = find(includechannels == header.ref_chan(1,i));
                if ~isempty(tmp)
                    ref_chan = [ref_chan tmp]; %#ok<AGROW>
                end
            end
        end
        header.ref_chan = ref_chan;
        clearvars ref_chan
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
if isfield(header,'channels')
    header.channels = header.channels(includechannels,:);
end

% Reduce data
data = data(includechannels,:,:);

% Correct header
if size(data,2) > 1
    header.numtimeframes = size(data,2);
end
if isfield(header,'includechans')
    header.includechans = header.includechans(1,includechannels);
else
    header.includechans = includechannels(:)';
end
if max(includechannels) <= header.numdatachannels | header.numdatachannels == header.numchannels
    header.numchannels = size(data,1);
    header.numdatachannels = header.numchannels;
    header.numauxchannels = 0;
    header.includedatachans = header.includechans;
else
    header.numchannels = size(data,1);
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
    if isfield(header.locs,'grad')
        header.locs.grad.label = header.locs.grad.label(includedatachannels,1);
        header.locs.grad.tra = header.locs.grad.tra(includedatachannels,:);
    end
    if isfield(header.locs,'sph_radius')
        header.locs.sph_radius = header.locs.sph_radius(1,includedatachannels);
        header.locs.sph_theta = header.locs.sph_theta(1,includedatachannels);
        header.locs.sph_phi = header.locs.sph_phi(1,includedatachannels);
    end
end

if isfield(header,'cov') & ~isempty(header.cov)
    includedatachannels = includechannels(1,1:end-header.numauxchannels);
    header.cov = header.cov(includedatachannels,includedatachannels);
end

return