% Perform dipol fit with results of ICA
%
% cfg = lab_ICA_dipolfit(activations,headeract,cfg)
%
% Written by F. Hatz 2013

function cfg = lab_ICA_dipolfit(activations,headeract,cfg)

if ~isfield(headeract,'W')
    return
end

if ~exist('cfg','var') | ~isfield(cfg,'Output_filepath')
    cfg.Output_file = headeract.EEG_file;
    cfg.Output_filepath = headeract.EEG_filepath;
end

if ~isfield(cfg,'ICABACK') | ~isfield(cfg.ICABACK,'IS') | isempty(cfg.ICABACK.IS)
    cfg.ICABACK = lab_set_dipolfit;
end

if ~exist('lab_inversesolution') %#ok<EXIST>
    disp('   Abort dipolfit, not implemented yet')
end

warning off %#ok<WNOFF>
mkdir(fullfile(cfg.Output_filepath,'DIPOLFIT'))
cfg.ICABACK.Output_filepath = fullfile(cfg.Output_filepath,'DIPOLFIT');
cfg.ICABACK.Output_file = cfg.Output_file;
warning on %#ok<WNON>

W = headeract.W;
for i = 1:size(activations,1)
    headeract.W = zeros(size(W));
    headeract.W(:,i) = W(:,i);
    [data,header] = lab_ica2eeg(activations,headeract);
    [RIStmp,cfg.ICABACK] = lab_inversesolution(data,header,cfg.ICABACK,true);
    if ~exist('RIS','var')
        RIS = RIStmp;
    else
        RIS.x(i,:) = RIS.x;
        RIS.y(i,:) = RIS.y;
        RIS.z(i,:) = RIS.z;
        RIS.data(i,:) = RIS.data;
        RIS.locs.x(1,i) = RIS.locs.x;
        RIS.locs.y(1,i) = RIS.locs.y;
        RIS.locs.z(1,i) = RIS.locs.z;
        RIS.locs.mom(:,i) = RIS.locs.mom;
    end
end
RIS.numsolutionpoints=size(RIS.data,1);
RIS.numchannels = size(RIS.data,1);
RIS.numdatachannels = size(RIS.data,1);
RIS.channels = char(num2str((1:size(data,1))'));
RIS.locs.labels = cellstr(num2str((1:size(RIS.data,1))'))';

for j = 1:size(cfg.ICABACK.IS.eformat,1)
    RIS_file = [cfg.Output_fileS '.' cfg.ICABACK.IS.eformat{j,1}];
    filename = fullfile(cfg.ICABACK.Output_filepath,RIS_file);
    if strcmp(cfg.ICABACK.IS.eformat{j,1},'edf') & size(RIS.data,1) > 200
        skipwriting = 1;
        disp(['   skip writing of ' cfg.ICABACK.IS.eformat{j,1} ' (to many channels)']);
    elseif ~strcmp(cfg.ICABACK.IS.eformat{j,1},'ris') & size(RIS.data,1) > 400
        skipwriting = 1;
        disp(['   skip writing of ' cfg.ICABACK.IS.eformat{j,1} ' (to many channels)']);
    else
        skipwriting = 0;
    end
    if skipwriting == 0
        if strcmp(cfg.ICABACK.IS.eformat{j,1},'ris')
            [~,cfg.IS.ICABACK] = lab_save_data(RIS,[],cfg.ICABACK.IS.eformat{j,1},filename,cfg.IS.ICABACK,cfg);
        else
            [~,cfg.IS.ICABACK] = lab_save_data(RIS.data,RIS,cfg.ICABACK.IS.eformat{j,1},filename,cfg.IS.ICABACK,cfg);
        end
    end
end