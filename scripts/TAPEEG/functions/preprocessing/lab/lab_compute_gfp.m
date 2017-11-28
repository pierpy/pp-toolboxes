% Calculate global field power and add the result as additional channel
%
% [data,header,gfp]=lab_compute_gfp(data,header)
%
% written by F. Hatz 2012

function [data,header,gfp] = lab_compute_gfp(data,header)

if isempty(data)
    return
end

if ~exist('header','var')
    header.numchannels = size(data,1);
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if ~isfield(header,'goodchans')
    header.goodchans = 1:header.numdatachannels;
end

if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
    disp('     Calculate gfp for channel reference')
    if isfield(header,'goodchans') & ~isempty(header.goodchans)
        avgref = mean(data(header.goodchans,:),1);
    else
        avgref = mean(data,1);
    end
    gfp = sqrt(sum((data - repmat(avgref,[size(data,1) 1])).^2,1)/size(data,1));
elseif isfield(header,'ref_chan') & ischar(header.ref_chan) & ...
        (strcmp(header.ref_chan,'mean') | strcmp(header.ref_chan,'median') | strcmp(header.ref_chan,'none'))
    disp(['     Calculate gfp for ' header.ref_chan  ' reference'])
    gfp = sqrt(sum(data.^2,1)/size(data,1));  
elseif ~isfield(header,'ref_chan')
    disp('     Calculate gfp without kown reference')
    gfp = sqrt(sum(data.^2,1)/size(data,1));  
else
    disp('     Abort: calculation of gfp not possible')
    gfp = [];
    return
end

data = cat(1,data,gfp);
header.numchannels = header.numchannels + 1;
if isfield(header,'numauxchannels')
    header.numauxchannels = header.numauxchannels + 1;
else
    header.numauxchannels = header.numchannels - header.numdatachannels;
end
if isfield(header,'channels')
    header.channels(end+1,1:3) = 'GFP';
end

end