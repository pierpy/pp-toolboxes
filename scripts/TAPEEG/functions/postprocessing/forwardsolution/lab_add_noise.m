% data = lab_add_noise(data,header,settings,novrb)
%
% data           = matrix (chans x numtimeframes)
% header         = output of lab_read_data
% settings.mode  = 1 (AR-model)
%                  2 (real data) -> dialog to open file
%                  3 (white noise)
%                  4 (pink noise)
% settings.coeff = coefficient for AR-model
% settings.dB    = dB of added noise
% novrb          = set to supress verbose-file
%
% written by F.Hatz 2013

function [data,header,cfg,skipprocessing] = lab_add_noise(data,header,cfg,novrb)

skipprocessing = 0;

if ~exist('novrb','var')
    novrb = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('novrb','var')
    if ~isfield(cfg,'EEG_file') & isfield(header,'EEG_filepath') & ~isempty(header.EEG_filepath)
        cfg.EEG_file = header.EEG_file;
        cfg.EEG_filepath = header.EEG_filepath;
        [~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
    else
        novrb = true;
    end
end

if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
    if ~isfield(cfg,'NOISE') | isempty(cfg.NOISE)
        [cfg,skipprocessing] = lab_set_add_noise(cfg,header);
        if skipprocessing == 1
            return
        end
    end
elseif ~isfield(cfg,'NOISE') | isempty(cfg.NOISE)
    return
end

if ~isfield(cfg.NOISE,'dB') | isempty(cfg.NOISE.dB)
    cfg.NOISE.dB = 10;
end
factor = (10^(cfg.NOISE.dB/10))^-1;

if isfield(header,'numdatachannels')
    NumChans = header.numdatachannels;
else
    NumChans = size(data,1);
end
Length = size(data,2);

if cfg.NOISE.mode == 1
    % add data from AR-model
    if ~isfield(cfg.NOISE,'coeff') | isempty(cfg.NOISE.coeff)
        cfg.NOISE.coeff = 10;
    end
    if exist('lab_eeg_ARmodel') %#ok<EXIST>
        noise = lab_eeg_ARmodel(NumChans,Length,cfg.NOISE.coeff);
    else
        disp('   Abort script ''lab_EEG_ARmodel'' missing')
    end
elseif cfg.NOISE.mode == 2
    % add real data
    if isfield(cfg.NOISE,'noise')
        noise = cfg.NOISE.noise;
    elseif isfield(cfg.NOISE,'data') & ~isempty(cfg.NOISE.data)
        noise = cfg.NOISE.data;
    elseif ~isfield(cfg,'MAIN') | ~isfiedl(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
        noise = lab_read_data;
    else
        noise = zeros(NumChans,Length);
    end
    if size(noise,1) >= Numchans
        noise = noise(1:Numchans,:);
    else
        noise = zeros(NumChans,Length);
    end
    if size(noise,2) >= Length
        noise = noise(:,1:Length);
    else
        disp('   Noise to short -> stitch several times')
        tmp = noise;
        headertmp = lab_create_header(tmp,true);
        headertmp.samplingrate = header.samplingrate;
        while size(noise,2) < Length
            noise = lab_stitch(noise,headertmp,tmp,headertmp,1);
        end
        noise = noise(:,1:Length);
    end
elseif cfg.NOISE.mode == 3
    % add white noise
    noise = 2*(rand(NumChans,Length)-0.5);
    noise = noise./max(abs(noise(:)));
elseif cfg.NOISE.mode == 4
    % add pink noise
    B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
    A = [1 -2.494956002 2.017265875 -0.522189400];
    noise = zeros(NumChans,Length);
    for i = 1:NumChans
        noise(i,:) = filter(B,A,2*(rand(1,Length)-0.5));
    end
    noise = noise./max(abs(noise(:)));
end

if ~exist('noise','var')
    disp('   Abort add noise, no valid mode selected')
    return
end

stddata = std(data(1:NumChans,:),[],2);
stdnoise = std(noise,[],2);
factor = (stddata ./ stdnoise) * factor;
for i = 1:size(noise,1)
    noise(i,:) = noise(i,:) * factor(i,1);
end
data(1:NumChans,:) = data(1:NumChans,:) + noise;

if novrb == false
    Verbose_file=[cfg.EEG_fileS '_addnoise.vrb'];
    fid=fopen(fullfile(cfg.EEG_filepath,Verbose_file),'w');
    fprintf(fid,'FAdd noise\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n');
    fprintf(fid,['Input file: ' cfg.EEG_file]);
    fprintf(fid,'\n\n');
    if cfg.NOISE.mode == 1
        fprintf(fid,'Mode: AR-model');
        fprintf(fid,'\n');
        fprintf(fid,['   Coefficient: ' num2str(cfg.NOISE.coeff)]);
    elseif cfg.NOISE.mode == 2
        fprintf(fid,'Mode: Real data');
    elseif cfg.NOISE.mode == 3
        fprintf(fid,'Mode: White noise');
    elseif cfg.NOISE.mode == 4
        fprintf(fid,'Mode: Pink noise');
    end
    fprintf(fid,'\n\n');
    fprintf(fid,['dB: ' num2str(cfg.NOISE.dB)]);
    fclose(fid);
end

end
