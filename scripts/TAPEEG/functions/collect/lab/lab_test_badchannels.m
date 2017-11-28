% Helper function for lab_collect_spectras
%
% [T,settings] = lab_test_badchannels(header,settings,T)
%
% written by F. Hatz 2012

function [T,settings] = lab_test_badchannels(header,settings,T)

if ~exist('settings','var')
    settings = [];
end

if ~isfield(settings,'QUALITY') | ~isfield(settings.QUALITY,'LAPL') | ~isfield(settings.QUALITY.LAPL,'lap_maxdistance')
    settings.QUALITY.LAPL.lap_maxdistance = 3;
end
settings.QUALITY.LAPL.lap_weightmaxdistance = 100;
if isfield(header,'locs')
    refmatrix = lab_calc_lapmatrix(header,settings.QUALITY.LAPL,1);
    refmatrix(refmatrix < 0) = 1;
elseif isfield(header,'numdatachannels')
    refmatrix = ones(header.numdatachannels,header.numdatachannels);
else
    refmatrix = ones(header.numchannels,header.numchannes);
end

T.badchans = zeros(1,size(refmatrix,1));
badtmp = zeros(1,size(refmatrix,1));

if isfield(settings.QUALITY,'badmode') & strcmp(settings.QUALITY.badmode,'interpolated') & ...
        isfield(header,'interpolated')
    if isfield(header,'interpolated') & ~isempty(header.interpolated)
        badtmp(header.interpolated) = 1;
    end
elseif isfield(header,'badchans') & ~isempty(header.badchans)
    badtmp(header.badchans) = 1;
end

for i = 1:size(refmatrix,1)
    T.badchans(1,i) = sum(badtmp(1,logical(refmatrix(i,:)))) / sum(refmatrix(i,:));
end

if isfield(header,'activationsexcluded')
    if isempty(header.activationsexcluded)
        T.badact(1,1) = 0;
    else
        T.badact(1,1) = mean(header.activationsexcluded);
    end
    T.badact(1,2) = length(badtmp) - sum(badtmp);
else
    T.badact = 0;
end