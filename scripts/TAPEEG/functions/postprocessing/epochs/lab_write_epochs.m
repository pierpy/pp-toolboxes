% Helper script for lab_save_epochs / Function to select epochs
%
% [epochs,cfg] = lab_write_epochs(epochs,header,cfg)
%
% Written by F. Hatz 2012

function [epochs,cfg] = lab_write_epochs(epochs,header,cfg)

if ~isfield(cfg,'EPOCH') | ~isfield(cfg.EPOCH,'folder')
    cfg.EPOCH.folder = 'Epochs';
end
if ~isfield(cfg.EPOCH,'eformat')
    cfg.EPOCH.eformat{1,1} = 'ep';
end
if ~isfield(cfg.EPOCH,'scaletxt')
    cfg.EPOCH.scaletxt = 'auto';
end
if ~isfield(cfg.EPOCH,'plotquality')
    cfg.EPOCH.plotquality = true;
end
if isfield(header,'events')
    header = rmfield(header,'events');
end

Epoch_filepath = fullfile(cfg.Output_filepath,cfg.EPOCH.folder);
warning off %#ok<WNOFF>
mkdir (fullfile(cfg.Output_filepath,cfg.EPOCH.folder));
warning on %#ok<WNON>

if cfg.EPOCH.minimalpart >= 1 & cfg.EPOCH.minimalpart <= size(epochs.data,3)
    maxepochs = cfg.EPOCH.minimalpart;
else
    maxepochs = size(epochs.data,3);
end

badoutsummary = [];
epochs.data = epochs.data(:,:,1:maxepochs);
epochs.goodmarker = epochs.goodmarkers(1,1:maxepochs);
epochs.badchans = epochs.badchans(1,1:maxepochs);
for i = 1:maxepochs
    filename = fullfile(Epoch_filepath,[cfg.Output_fileS '_E' num2str(i) '.sef']);
    headertmp = header;
    if isfield(cfg.EPOCH,'replacebad') & cfg.EPOCH.replacebad == true
        if isfield(epochs,'badchans')
            headertmp.badchans = epochs.badchans{1,i};
            headertmp.goodchans = setdiff(1:headertmp.numdatachannels,headertmp.badchans);
            headertmp.goodchans = headertmp.goodchans(:)';
        end
        if isfield(epochs,'interpolated')
            headertmp.interpolated = epochs.interpolated{1,i};
        end
    end
    headertmp.numtimeframes = size(epochs.data,2);
    if isfield(epochs,'goodsum') & length(epochs.goodsum) >= i
        headertmp.quality = epochs.goodsum(i);
    end
    for j = 1:length(cfg.EPOCH.eformat)
        [filename,cfg.EPOCH] = lab_save_data(epochs.data(:,:,i),headertmp,cfg.EPOCH.eformat{j},filename,cfg.EPOCH,cfg);
    end
    badouttitle{i,1} = [cfg.Output_fileS '_E' num2str(i)]; %#ok<AGROW>
    tmp = {badouttitle{i,1},num2str(length(epochs.badchans{1,i}));'Spectral50Hz',num2str(epochs.bad.Spect50{1,i}); ...
        'Spectral60Hz',num2str(epochs.bad.Spect60{1,i});'SpectralLow',num2str(epochs.bad.SpectLow{1,i}); ...
        'SpectralHigh',num2str(epochs.bad.SpectHigh{1,i});'Variance',num2str(epochs.bad.Var{1,i}); ...
        'Hurst',num2str(epochs.bad.Hurst{1,i});'Marker',num2str(epochs.bad.Marker{1,i});'',''};
    badoutsummary = cat(1,badoutsummary,tmp);
    clearvars tmp
end

% Write Bad chans info
badoutheader{1,1} = 'file';
for i = 1:size(epochs.data,1)
    badoutheader{1,i+1} = header.channels(i,:);
end
badout = cat(2,badouttitle,num2cell(epochs.badoutelectrodes(1:size(badouttitle,1),:)));
badout = cat(1,badoutheader(1,1:size(badout,2)),badout);
dirtmp = cd(Epoch_filepath);
lab_write_xls([cfg.Output_fileS '_Bad.xls'],badout');
lab_write_xls([cfg.Output_fileS '_BadInfo.xls'],badoutsummary);
if isfield(header,'badchans') & ~isempty(header.badchans) & header.badchans ~= 0
    badout2 = zeros(1,header.numchannels);
    badout2(1,header.badchans) = 1;
    badout2 = num2cell(badout2);
    badout2 = ['all epochs' badout2];
    badout2 = cat(1,badoutheader,badout2);
    lab_write_xls([cfg.Output_fileS '_Bad_IP.xls'],badout2');
end

% Write Marker-file
if isfield(epochs,'events')
    Marker_file=[cfg.Output_fileS '.mrk'];
    lab_write_mrk(Marker_file,epochs);
end

% Save quality check plots
if isfield(header,'locs') & cfg.EPOCH.plotquality == true
    lab_epoch_qualitycontrol(epochs,header,fullfile(Epoch_filepath,cfg.Output_fileS));
end

cd(dirtmp);

return
