% Interpolate bad channels
%
% [data,header,method] = lab_interpolate_bad(data,header,method)
%
% data   = matrix (chans x timeframes)
% header = output of lab_read_data
% method = '3D' (3D interpolation with TriScatterdInterp)
%          'spherical' (spherical spline interpolation - code from eeglab)
%
% written by F. Hatz 2012

function [data,header,method] = lab_interpolate_bad(data,header,method,noverbose)

if ~exist('method','var') | isempty(method)
    method = 'spherical';
end
if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end

if size(data,3) > 1
    nepochs = size(data,3);
    data = reshape(data,[size(data,1) size(data,2)*size(data,3)]);
end

if isfield(header,'badchans') & ~isempty(header.badchans) & ~header.badchans == 0 & isfield(header,'locs') & ~isfield(header.locs,'grad')
    if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
        header.badchans = setdiff(header.badchans,header.ref_chan);
        header.goodchans = setdiff(1:size(header.locs.x,2),union(header.ref_chan,header.badchans));
    else
        header.goodchans = setdiff(1:size(header.locs.x,2),header.badchans);
    end
    header.badchans = header.badchans(:)';
    header.goodchans = header.goodchans(:)';
    if ~isfield(header,'splitchans')
        if strcmp(method,'3D')
            if ~exist('noverbose','var')
                disp('     3-D interpolation bad channels');
            end
            [xbad ,ybad] = pol2cart(header.locs.theta(header.badchans),header.locs.radius(header.badchans));
            [xgood ,ygood] = pol2cart(header.locs.theta(header.goodchans),header.locs.radius(header.goodchans));
            badtmp = zeros(length(header.badchans),size(data,2));
            for t = 1:size(data,2)
               F = TriScatteredInterp(ygood', xgood', data(header.goodchans,t),'nearest');
                badtmp(:,t) = F(ybad', xbad');
                clearvars F;
            end
            data(header.badchans,:) = badtmp;
            clearvars badtmp xbad ybad xgood ygood;
            if isfield(header,'interpolated')
                header.interpolated = union(header.interpolated,header.badchans);
            else
                header.interpolated = header.badchans;
            end
            header.goodchans = union(header.goodchans,header.badchans);
            header.badchans = [];
        else
            if ~exist('noverbose','var')
                disp('     spherical-spline-Interpolation bad channels')
            end
            xelec = header.locs.x(header.goodchans);
            yelec = header.locs.y(header.goodchans);
            zelec = header.locs.z(header.goodchans);
            rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
            xelec = xelec./rad;
            yelec = yelec./rad;
            zelec = zelec./rad;
            xbad = header.locs.x(header.badchans);
            ybad = header.locs.y(header.badchans);
            zbad = header.locs.z(header.badchans);
            rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
            xbad = xbad./rad;
            ybad = ybad./rad;
            zbad = zbad./rad;
            [~ , ~ , ~, badchansdata] = lab_spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, data(header.goodchans,:));
            data(header.badchans,:) = badchansdata;
            header.interpolated = header.badchans;
            header.goodchans = union(header.goodchans,header.badchans);
            header.badchans = [];
            clearvars xelec yelec zelec xbad ybad zbad rad badchansdata P;
        end
    else
        if strcmp(method,'3D')
            if ~exist('noverbose','var')
                disp('     3-D interpolation bad channels');
            end
            [x,y]  = pol2cart(header.locs.theta,header.locs.radius);
            for j = 1:size(header.splitchans,1)
                badchans = intersect(header.goodchans,(header.splitchans(j,1):header.splitchans(j,2)));
                goodchans = intersect(header.goodchans,(header.splitchans(j,1):header.splitchans(j,2)));
                for t = 1:size(data,2)
                    F = TriScatteredInterp(y(goodchans)', x(goodchans)', data(goodchans,t),'nearest');
                    data(badchans,:) = F(y(badchans)', x(badchans)');
                    clearvars F;
                end
                clearvars goodchans badchans
            end
            clearvars badtmp x y
            header.interpolated = header.badchans;
            header.goodchans = union(header.goodchans,header.badchans);
            header.badchans = [];
        else
            if ~exist('noverbose','var')
                disp('     spherical-spline-Interpolation bad channels')
            end
            xelec = header.locs.x;
            yelec = header.locs.y;
            zelec = header.locs.z;
            rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
            xelec = xelec./rad;
            yelec = yelec./rad;
            zelec = zelec./rad;
            xbad = header.locs.x(header.badchans);
            ybad = header.locs.y(header.badchans);
            zbad = header.locs.z(header.badchans);
            rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
            xbad = xbad./rad;
            ybad = ybad./rad;
            zbad = zbad./rad;
            for i = 1:size(xbad,2)
                cluster = [];
                for j = 1:size(header.splitchans,1)
                    if find((header.splitchans(j,1):header.splitchans(j,2)) == header.badchans(i))
                        cluster = j;
                    end
                end
                if ~isempty(cluster)
                    goodchans = intersect(header.goodchans,(header.splitchans(cluster,1):header.splitchans(cluster,2)));
                    if ~isempty(goodchans)
                        [~ , ~ , ~, badchansdata(i,:)] = lab_spheric_spline(xelec(goodchans),yelec(goodchans), ...
                            zelec(goodchans),xbad(1,i),ybad(1,i),zbad(1,i),data(goodchans,:)); %#ok<AGROW>
                    else
                        badchansdata(i,:) = zeros(1,size(data,2)); %#ok<AGROW>
                    end
                    clearvars goodchans
                else
                    badchansdata(i,:) = zeros(1,size(data,2)); %#ok<AGROW>
                end
            end
            data(header.badchans,:) = badchansdata;
            if isfield(header,'interpolated')
                header.interpolated = union(header.interpolated,header.badchans);
            else
                header.interpolated = header.badchans;
            end
            header.goodchans = union(header.goodchans,header.badchans);
            header.badchans = [];
            clearvars xelec yelec zelec xbad ybad zbad rad badchansdata P cluster;
        end
    end
    if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
        header.badchans = setdiff(header.badchans,header.ref_chan);
    end
    header.goodchans = setdiff(1:header.numdatachannels,header.badchans);
    header.badchans = header.badchans(:)';
    header.goodchans = header.goodchans(:)';
elseif ~isfield(header,'locs') | isfield(header.locs,'grad')
    disp('   Abort: Interpolation of bad channels not possible, loc-file missing')
end

if exist('nepochs','var')
    data = reshape(data,[size(data,1) size(data,2)/nepochs nepochs]);
    header.numtimeframes = size(data,2);
end
