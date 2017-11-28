% Recalculate EEG with result of ICA
%
% [data,header] = lab_ica2eeg(activations,header)
%
% written by F. Hatz 2012

function [data,header] = lab_ica2eeg(activations,header)

if ~isfield(header,'W')
    disp('   Abort: no ICA information')
else
    W = header.W;
    header = rmfield(header,'W');
end
if isfield(header,'ICAchans')
    ICAchans = header.ICAchans;
    header = rmfield(header,'ICAchans');
else
    ICAchans = 1:size(activations,1);
end
if isfield(header,'ICAmaxchan')
    data = zeros(header.ICAmaxchan,size(activations,2));
    header = rmfield(header,'ICAmaxchan');
else
    data = zeros(max(ICAchans),size(activations,2));
end
    
if isfield(header,'badchans') & ~isempty(header.badchans) & min(header.badchans) > 0
    W(:,header.badchans)=0;
end
data(ICAchans,:)=W*activations;
header.activationsexcluded = length(activations_exclude);
header.numchannels = size(data,1);
header.numdatachannels = size(data,1);
if isfield(header,'ICAref_chan')
    header.ref_chan = header.ICAref_chan;
    header = rmfield(header,'ICAref_chan');
    if isnumeric(header.ref_chan)
        ICAchans = union(ICAchans,header.ref_chan);
    end
else
    header.ref_chan = 'unkown';
end
header.goodchans = ICAchans;
if size(header.goodchans,1) > 1
    header.goodchans = header.goodchans';
end
header.badchans = setdiff(1:size(data,1),ICAchans);
header.badchans = header.badchans(:)';

end