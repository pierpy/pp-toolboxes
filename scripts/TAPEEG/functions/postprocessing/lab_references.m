% Calculate new references
%
% [data,header,cfg] = lab_references(data,header,method,settings)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% method   = mean (only good channels included)
%            median (only good channels included)
%            laplacian (= reference for every channel is gaussian
%                       distribution of surrounding electrodes)
%            montage(= defined by xls-file with montage information)
%            channels (single or multiple channel numbers)
% settings = structure with config (optional)
%
% written by F. Hatz 2012

function [data,header,settings] = lab_references(data,header,method,settings)

if ~exist('settings','var')
    settings = [];
end
if ~exist('header','var') | ~exist('method','var')
    disp('   Re-referencing not possible (header or method is missing)');
    return
end

if ischar(method) & ~isnan(str2double(method))
    method = str2num(method); %#ok<ST2NM>
end
if isstruct(method)
    if length(method) > 1
        disp(['   Process only first montage: ' method(1).name]);
        method = method(1);
    end
    disp(['   Compute Montage: ' method.name]);
elseif isnumeric(method)
    disp(['   Compute ' num2str(method) '-reference']);
elseif ischar(method)
    disp(['   Compute ' method '-reference']);
else
    disp('   No valid method for re-reference');
    return
end

if size(data,3) > 1
    nepochs = size(data,3);
    data = reshape(data,[size(data,1) size(data,2)*size(data,3)]);
end
if ~(header.numtimeframes == size(data,2))
    header.numtimeframes = size(data,2);
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end
if ~isfield(header,'ref_chan')
    header.ref_chan = [];
elseif isnumeric(header.ref_chan) & max(header.ref_chan) > header.numdatachannels
    disp('   omit reference channel as index of reference channel is higher than number of channels')
    header.ref_chan = [];
end
if ~isfield(header,'goodchans')
    header.badchans = [];
    header.goodchans = 1:header.numdatachannels;
end
if header.badchans == 0
    header.badchans = [];
end
if isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & min(header.ref_chan) > 0
    header.badchans = setdiff(header.badchans,header.ref_chan);
    header.badchans = header.badchans(:)';
    header.goodchans = setdiff(header.goodchans,header.ref_chan);
    header.goodchans = header.goodchans(:)';
end
if ischar(header.ref_chan) & ~isempty(header.ref_chan)
    refchar = true;
else
    refchar = false;
end

%--------------------------------------------------------------------------
% Compute EEG against mean reference
%--------------------------------------------------------------------------
if strcmp(method,'mean')
    if refchar == true
        disp(['   Reference is ' header.ref_chan ' calculation of mean not possible']);
    else
        if isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & min(header.ref_chan) > 0
            goodchans = setdiff(header.goodchans,header.ref_chan);
            header.badchans = setdiff(header.badchans,header.ref_chan);
            header.badchans = header.badchans(:)';
            header.goodchans = union(goodchans,header.ref_chan);
            header.goodchans = header.goodchans(:)';
        else
            goodchans = header.goodchans;
        end
        if max(header.ref_chan) > header.numdatachannels
            header.ref_chan = [];
        end
        refmatrix = zeros(header.numdatachannels,header.numdatachannels);
        refmatrix(:,goodchans) = -1 / length(goodchans);
        refmatrix(1:header.numdatachannels+1:end) = refmatrix(1:header.numdatachannels+1:end) + 1;
        % -- old version with uncorrect mean value --
        % badchans = setdiff(1:header.numdatachannels,goodchans);
        % refmatrix(goodchans,goodchans) = -1 / (length(goodchans) - 1);
        % refmatrix(badchans,goodchans) = -1 / (length(goodchans));
        % refmatrix(1:header.numdatachannels+1:end) = 1;
        % -------------------------------------------
        data(1:header.numdatachannels,:) = refmatrix * data(1:header.numdatachannels,:);
        header.ref_chan = 'mean';
        clearvars refmatrix
    end
end

%--------------------------------------------------------------------------
% Compute EEG against median reference
%--------------------------------------------------------------------------
if strcmp(method,'median')
    if refchar == true
        disp(['   Reference is ' header.ref_chan ' calculation of median not possible']);
    else
        if isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & min(header.ref_chan) > 0
            goodchans = setdiff(header.goodchans,header.ref_chan);
            header.badchans = setdiff(header.badchans,header.ref_chan);
            header.badchans = header.badchans(:)';
            header.goodchans = union(goodchans,header.ref_chan);
            header.goodchans = header.goodchans(:)';
        else
            goodchans = header.goodchans;
        end
        reference = zeros(header.numdatachannels,size(data,2));
        for i = 1:header.numdatachannels
            reference(i,:) = median(data(setdiff(goodchans,i),:),1);
        end
        data(1:header.numdatachannels,:) = data(1:header.numdatachannels,:) - reference;
        header.ref_chan = 'median';
        clearvars reference
    end
end

%--------------------------------------------------------------------------
% Compute EEG against Laplacian reference
%--------------------------------------------------------------------------
if strcmp(method,'laplacian')
    if isfield(header,'locs')
        if ~exist('settings','var') | ~isfield(settings,'LAPL')
            settings.LAPL.lap_maxdistance=4;
            settings.LAPL.lap_weightmaxdistance=5;
        end
        [refmatrix,sigma,settings.LAPL] = lab_calc_lapmatrix(header,settings.LAPL);
        disp('     compute laplacian-referenced data');
        data(1:header.numdatachannels,:) = refmatrix * data(1:header.numdatachannels,:);
        header.ref_lapl.sigma = sigma;
        header.ref_lapl.maxdistance = settings.LAPL.lap_maxdistance;
        header.ref_lapl.weightmaxdistance = settings.LAPL.lap_weightmaxdistance;
        if isfield(settings.LAPL,'lap_excluderef') & settings.LAPL.lap_excluderef == true & ...
                isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & min(header.ref_chan) > 0
            [data,header] = lab_reduce_channels(data,header,setdiff(1:size(data,1),header.ref_chan));
        elseif isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & min(header.ref_chan) > 0
            header.goodchans = union(header.goodchans,header.ref_chan);
            header.goodchans = header.goodchans(:)';
        end
        header.ref_chan = 'laplacian';
        clearvars distances refmatrix i j sigma
    else
        disp('     calculating laplacian-reference not possible (missing locs)')
    end
end
%--------------------------------------------------------------------------
% Compute EEG against selected electrodes
%--------------------------------------------------------------------------
if isnumeric(method) & max(method) <= header.numdatachannels
    if length(method) > 1
        reference = mean(data(method,:),1);
    else
        reference = data(method,:);
    end
    for i = 1:header.numdatachannels
        data(i,:) = data(i,:) - reference;
    end
    if isnumeric(header.ref_chan) & ~isempty(header.ref_chan) & min(header.ref_chan) > 0
        header.goodchans = union(header.goodchans,header.ref_chan);
        header.goodchans = header.goodchans(:)';
    end
    header.ref_chan = method;
    clearvars reference
end
%--------------------------------------------------------------------------
% Compute Montage
%--------------------------------------------------------------------------
if isstruct(method) & isfield(method,'chans')
    if isfield(method,'numchans') & method.numchans > header.numchannels
        disp('    No referencing, montage information not matching')
        return
    end
    datatmp = zeros(size(method.chans,1),size(data,2));
    flag = false(1,size(method.chans,1));
    for i = 1:size(method.chans,1);
        if ~isempty(method.chans{i,1}(1,1))
            flag(i) = true;
            active = method.chans{i,1}(1,1) + (method.chans{i,2} * header.numdatachannels);
            if isnumeric(method.chans{i,3}) | isempty(method.chans{i,3})
                if isempty(method.chans{i,3}) | isnan(method.chans{i,3})
                    reference = 0;
                else
                    reference = method.chans{i,3} + (method.chans{i,4} * header.numdatachannels);
                end
                if active <= header.numchannels & max(reference) <= header.numchannels & active > 0
                    if reference == 0
                        datatmp(i,:) = data(active,:);
                    else
                        if size(reference,2) > 1
                            datatmp(i,:) = data(active,:) - mean(data(reference,:),1);
                        else
                            datatmp(i,:) = data(active,:) - data(reference,:);
                        end
                    end
                end
            elseif strcmp(method.chans{i,3},'AVG')
                datatmp(i,:) = data(active,:) - mean(data(header.goodchans,:),1);
            elseif strcmp(method.chans{i,3},'LAPL')
                if active <= header.numdatachannels
                    if ~exist('settings','var')
                        settings = [];
                    end
                    [refmatrix,sigma] = lab_calc_lapmatrix(header,settings);
                    datatmp(i,:) = refmatrix(active,:) * data(1:header.numdatachannels,:);
                else
                    datatmp(i,:) = data(active,:);
                end
            end
        end
    end
    if max(flag) == 0
        data = [];
        header = [];
        return
    end
    data = datatmp;
    header.numchannels = size(data,1);
    header.numdatachannels = size(data,1);
    header.ref_chan = method.name;
    header.channels = char(method.label);
    header.numauxchannels = 0;
    header.ecg_ch = 0;
    if isfield(header,'locs')
        header = rmfield(header,'locs');
    end
    if exist('sigma','var')
        header.ref_lapl.sigma = sigma;
        header.ref_lapl.maxdistance = settings.LAPL.lap_maxdistance;
        header.ref_lapl.weightmaxdistance = settings.LAPL.lap_weightmaxdistance;
    end
    header.montage = method;
    header.montage = true;
    if isfield(header,'interpolated')
        header = rmfield(header,'interpolated');
    end
    header.badchans = [];
    header.goodchans = 1:header.numchannels;
end

if exist('nepochs','var')
    data = reshape(data,[size(data,1) size(data,2)/nepochs nepochs]);
    header.numtimeframes = size(data,2);
end

return

