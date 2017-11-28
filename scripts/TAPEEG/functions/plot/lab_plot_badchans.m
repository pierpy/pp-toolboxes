% Plot electrodes with highlighted bad channels
%
% lab_plot_badchans(header)
%
%      header     see lab_create_header
%
% written by F. Hatz 2013

function badchans = lab_plot_badchans(header)

if ~isfield(header,'locs') | isempty(header.locs)
    disp('no LOCS')
    badchans = [];
    return
end
if ~isfield(header,'badchans') | isempty(header.badchans)
    disp('no bad channels')
    badchans = [];
    return
end

settings.LOCS = header.locs;
settings.indexed = header.badchans;
settings.Color = [1 1 1];
settings.ColorIdx = [1 0 0];

badchans = lab_plot_locs(settings,1,0,0);