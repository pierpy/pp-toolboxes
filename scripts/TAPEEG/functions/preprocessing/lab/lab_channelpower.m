% Calculate channel power
%
% [data,header,power]=lab_channelpower(data,header)
%
% written by F. Hatz 2014

function [data,header,power]=lab_channelpower(data,header)

if isempty(data)
    return
end

if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end

power = sqrt(sum(data(1:header.numdatachannels,:).^2,2) / size(data(1:header.numdatachannels),2));