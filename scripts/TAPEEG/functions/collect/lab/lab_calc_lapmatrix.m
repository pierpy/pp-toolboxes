% Calculate laplacian distance-matrix for electrode-channels
%
% [refmatrix,sigma,LAPL] = lab_calc_lapmatrix(header,LAPL,novrb)
%
%    header                      see lab_create_header(data)
%    LAPL.lap_maxdistance        maximal distance in minimal distances between
%                                electrodes (default = 4)
%    LAPL.lap_weightmaxdistance  weight of electrodes at maximal distance
%                                in percent (default = 5)
%    novrb                       1 = silent mode
%
% written by F. Hatz 2012

function [refmatrix,sigma,LAPL] = lab_calc_lapmatrix(header,LAPL,novrb)

global Settings_Path

if ~exist('novrb','var')
    novrb = 0;
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end

% Read electrodes
if ~isfield(header,'locs')
    if novrb == 0
        disp ('   Read electrodes position file')
    end
    if exist(fullfile(Settings_Path,'electrodes.sfp'),'file')
        electrodesfile = fullfile(Settings_Path,'electrodes.sfp');
    elseif exist(fullfile(Settings_Path,'electrodes.els'),'file')
        electrodesfile = fullfile(Settings_Path,'electrodes.els');
    elseif exist(fullfile(Settings_Path,'electrodes.xyz'),'file')
        electrodesfile = fullfile(Settings_Path,'electrodes.xyz');
    end
    if exist('electrodesfile','var')
        header.locs = lab_read_locs(electrodesfile);
        if isfield(header,'includedatachans')
            header.locs.x = header.locs.x(1,header.includedatachans);
            header.locs.y = header.locs.y(1,header.includedatachans);
            header.locs.z = header.locs.z(1,header.includedatachans);
            header.locs.theta = header.locs.theta(1,header.includedatachans);
            header.locs.radius = header.locs.radius(1,header.includedatachans);
        end
    end
end

if isfield(header,'locs') & length(header.locs.x) == header.numdatachannels
    % Calculate distance matrix
    if novrb == 0
        disp ('     calculate distance matrix')
    end
    distances = zeros(length(header.locs.x),length(header.locs.x));
    for i = 1:length(header.locs.x);
        for j = 1:length(header.locs.x);
            distances(i,j) = ((header.locs.x(1,i) - header.locs.x(1,j))^2 + ...
                (header.locs.y(1,i) - header.locs.y(1,j))^2 + ...
                (header.locs.z(1,i) - header.locs.z(1,j))^2)^0.5;
        end
    end
    clearvars i j
    if novrb == 0
        disp('     normalize distance (minimal distance = 1)')
    end
    distances = distances / min(distances(distances > 0));
    
    %Calculate refmatrix weights
    if ~exist('LAPL','var') | ~isfield(LAPL,'lap_maxdistance')
        LAPL.lap_maxdistance = 4;
        LAPL.lap_weightmaxdistance = 5;
    end
    refmatrix = zeros(length(header.locs.x),length(header.locs.x));
    sigma =   (LAPL.lap_maxdistance^2 / (-2 * log(LAPL.lap_weightmaxdistance / 100)))^0.5;
    for i = 1:length(header.locs.x);
        for j = 1:length(header.locs.x);
            if distances(i,j) <= LAPL.lap_maxdistance
                if ~isinf(sigma)
                    refmatrix(i,j) = exp(-distances(i,j)^2 / (sigma^2 * 2));
                else
                    refmatrix(i,j) = 1;
                end
            else
                refmatrix(i,j) = 0;
            end
        end
    end
    if isfield(header,'badchans') & ~isempty(header.badchans)
        refmatrix(:,header.badchans) = 0;
    end
    if isnumeric(header.ref_chan) & ~isempty(header.ref_chan)
        refmatrix(:,header.ref_chan) = 0;
    end
    refmatrix(1:size(refmatrix,1)+1:end) = 0;
    refmatrix =   - refmatrix ./ repmat(( sum(refmatrix,2)),1,size(refmatrix,1));
    refmatrix(1:size(refmatrix,1)+1:end) = 1;
    refmatrix(isnan(refmatrix)) = 0;
else
    disp('     error: laplacian reference not possible, electrodes locations unkown')
    refmatrix = zeros(header.numdatachannels,header.numdatachannels);
    refmatrix(1:size(refmatrix,1)+1:end) = 1;
end

return