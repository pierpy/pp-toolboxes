% (ICA) independend component analysis
%
% [activations,W,cfg,badAll] = lab_ICAstart(data,header,cfg)
%
% data     = matrix (chans x timeframes)
% header   = output of lab_read_data
% cfg      = structure with config (optional)
%
% written by F. Hatz 2012

function [activations,headerICA,cfg] = lab_ICAstart(data,header,cfg)

disp('ICA')

if ~exist('header','var') | ~isfield(header,'samplingrate')
    header = lab_create_header(data);
end
if ~exist('cfg','var')
    cfg = [];
end

if~isfield(cfg,'EEG_file')
    cfg.EEG_file = header.EEG_file;
    cfg.EEG_filepath = header.EEG_filepath;
end
[~,~,~,cfg.EEG_fileS] = lab_filename(cfg.EEG_file);
EEG_file = cfg.EEG_fileS;

if ~isfield(header,'goodchans') | isempty(header.goodchans)
    header.goodchans = 1:header.numdatachannels;
    header.badchans = [];
end

if ~isfield(cfg,'ICA') | ~isfield(cfg.ICA,'good')
    [cfg,skipprocessing] = lab_set_ICAstart(cfg,header,false);
    pause(0.2);
    if skipprocessing == 1
        activations = [];
        headerICA = [];
        return
    end
end


%--------------------------------------------------------------------------
% reduce data length if necessary
%--------------------------------------------------------------------------
if isfield(cfg,'ICA') & isfield(cfg.ICA,'epoch') & cfg.ICA.epoch > 0
    if length(cfg.ICA.epoch) == 1
        if cfg.ICA.epoch + 5*header.samplingrate <= size(data,2)
            cfg.ICA.epoch = [5*header.samplingrate (cfg.ICAepoch+5*header.samplingrate)];
        elseif cfg.ICA.epoch <= size(data,2)
            cfg.ICA.epoch = [1 cfg.ICA.epoch];
        else
            cfg.ICA.epoch = [1 size(data,2)];
        end
    elseif length(cfg.ICA.epoch) > 2 & cfg.ICA.epoch(1) < size(data,2) & cfg.ICA.epoch(2) <= size(data,2)
        cfg.ICA.epoch = cfg.ICA.epoch(1,1:2);
    else
        cfg.ICA.epoch = [1 size(data,2)];
    end
    settings.REDUCE.firstsample = cfg.ICA.epoch(1);
    settings.REDUCE.lastsample = cfg.ICA.epoch(2);
    [data,header] = lab_reduce_datalength(data,header,settings);
end

%--------------------------------------------------------------------------
% define channels to include in ICA
%--------------------------------------------------------------------------
if cfg.ICA.good == 1;
    if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
        header.goodchans = setdiff(header.goodchans,cfg.exclude);
        header.goodchans = header.goodchans(:)';
    end
    ICAchans = header.goodchans;
else
    ICAchans = 1:header.numdatachannels;
    if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
        ICAchans = setdiff(ICAchans,header.ref_chan);
    end
    if isfield(cfg,'exclude') & ~isempty(cfg.exclude)
        ICAchans = setdiff(ICAchans,cfg.exclude);
    end
end

% Look for channels without data
eeg_chans = 1:header.numdatachannels;
for ch=1:header.numdatachannels
    tmp(ch)=std(data(ch,:)); %#ok<AGROW>
end
chanNoSignal = tmp<10^-14;
chanNoSignal = eeg_chans(logical(chanNoSignal));
if min(chanNoSignal) > 0
    ICAchans = setdiff(ICAchans,chanNoSignal);
end
clearvars tmp eeg_chans chanNoSignal

%--------------------------------------------------------------------------
% Create Output-directory
%--------------------------------------------------------------------------
warning off %#ok<WNOFF>
cfg.ICA_filepath = fullfile(cfg.EEG_filepath,'ICA');
mkdir (cfg.ICA_filepath);
warning on %#ok<WNON>
cfg.Output_filepath = cfg.ICA_filepath;
cfg.Output_file = [EEG_file '_ICA.sef'];
cfg.ICA_file=[EEG_file '_ICA.sef'];

%--------------------------------------------------------------------------
% Write EEG before ICA
%--------------------------------------------------------------------------
disp('   Write file before ICA')
cfg.EEG_output_file=[EEG_file 'filt.sef'];
lab_write_sef(fullfile(cfg.ICA_filepath,cfg.EEG_output_file),data,header);

%--------------------------------------------------------------------------
% ICA decomposition
%--------------------------------------------------------------------------
if ~isfield(cfg.ICA,'pca')
    cfg.ICA.pca = 0;
elseif cfg.ICA.pca == -1
    cfg.ICA.pca = rank(data(ICAchans,:));
end
if ~isfield(cfg.ICA,'extended')
    cfg.ICA.extended = 0;
end
if ~isfield(cfg.ICA,'type')
    cfg.ICA.type = 'runica';
end
disp('   Start ICA')
warning off %#ok<WNOFF>
W = [];
dorank = 0;
while isempty(W)
    switch cfg.ICA.type
        case 'runica'
            if cfg.ICA.pca > 0
                [weights, sphere, ~, ~, ~, ~, activations] = runica(data(ICAchans,:),'verbose','on','extended',cfg.ICA.extended,'pca',cfg.ICA.pca);
            else
                [weights, sphere, ~, ~, ~, ~, activations] = runica(data(ICAchans,:),'verbose','on','extended',cfg.ICA.extended);
            end
            W = pinv(weights*sphere);
        case 'binica'
            if cfg.ICA.pca > 0
                [weights, sphere] = runica(data(ICAchans,:),'verbose','on','extended',cfg.ICA.extended,'pca',cfg.ICA.pca);
            else
                [weights, sphere] = runica(data(ICAchans,:),'verbose','on','extended',cfg.ICA.extended);
            end
            activations = weights * sphere * data(ICAchans,:);
            W = pinv(weights*sphere);
        case 'fastICA'
            [activations,W] = fastica(data(ICAchans,:),'displayMode','off');
            clearvars tmp
        case 'jader'
            if cfg.ICA.pca > 0
                [weights] = jader(data(ICAchans,:),cfg.ICA.pca);
            else
                [weights] = jader(data(ICAchans,:));
            end
            sphere  = eye(size(weights,2));
            activations = weights * sphere * data(ICAchans,:);
            W = pinv(weights*sphere);
    end
    if ~isreal(W)
        if isempty(cfg.ICA.pca) | cfg.ICA.pca == 0
            if dorank == 0
                dorank = 1;
                cfg.ICA.pca = rank(data(ICAchans,:));
                if cfg.ICA.pca >= length(ICAchans)
                    cfg.ICA.pca = length(ICAchans) - 1;
                end
                W = [];
                disp(['     ICA error, reduce number of components to ' num2str(cfg.ICA.pca)])
            else
                cfg.ICA.pca = length(ICAchans) - 1;
                W = [];
                disp(['     ICA error, reduce number of components to ' num2str(cfg.ICA.pca)])
            end
        elseif cfg.ICA.pca > 1
            cfg.ICA.pca = cfg.ICA.pca - 1;
            W = [];
            disp(['     ICA error, reduce number of components to ' num2str(cfg.ICA.pca)])
        end
    end
end
warning on %#ok<WNON>

% normalize weigths of activations
disp('     normalize weigths of activations')
tmp = max(abs(W),[],1);
W = W ./ repmat(tmp,size(W,1),1);
activations = activations .* repmat(tmp',1,size(activations,2));

% sort activations by weight
disp('     sort activations by weight')
tmp = sum(abs(W),1);
[~,tmp] = sort(tmp,'descend');
W = W(:,tmp);
activations = activations(tmp,:);

topo=zeros(size(W,2),header.numdatachannels);
topo(:,ICAchans)=W';
topogood = W';
clearvars weights compvars signs bias lrates sphere;
headerICA = lab_reduce_header(header,ICAchans);
headerICA.numtimeframes = size(activations,2);
headerICA.numchannels = size(activations,1);
headerICA.numdatachannels = headerICA.numchannels;
headerICA.channels = num2str((1:size(activations,1))');
headerICA.numauxchannels = 0;
headerICA.W = W;
if isfield(header,'bad')
    headerICA = rmfield(headerICA,'bad');
end
headerICA.EEG_file = cfg.ICA_file;
headerICA.EEG_filepath = cfg.ICA_filepath;

%--------------------------------------------------------------------------
% Detect bad components
%--------------------------------------------------------------------------
if isfield(cfg.ICA,'BAD') & ~isempty(cfg.ICA.BAD)
    disp('   Detect bad activations')
    if isfield(cfg.ICA.BAD,'ecgdetect') & cfg.ICA.BAD.ecgdetect == 1 & cfg.ICA.BAD.ecg_ch <= size(data,1)
        cfg.ICA.BAD.ecg_chan = data(cfg.ICA.BAD.ecg_ch,:);
    end
    if isfield(cfg.ICA.BAD,'eog') & ~isempty(cfg.ICA.BAD.eog)
        cfg.ICA.BAD = lab_calculate_eog(data,header,cfg.ICA.BAD);
    end
    tEEG_file = cfg.EEG_file;
    tEEG_filepath = cfg.EEG_filepath;
    cfg.EEG_file = cfg.ICA_file;
    cfg.EEG_filepath = cfg.ICA_filepath;
    if isfield(headerICA,'locs')
        [badAll,bad,cfg.ICA.BAD] = lab_detect_bad(activations,rmfield(headerICA,'locs'),cfg.ICA.BAD,cfg);
    else
        [badAll,bad,cfg.ICA.BAD] = lab_detect_bad(activations,headerICA,cfg.ICA.BAD,cfg);
    end
    bad.All = badAll;
    cfg.EEG_file = tEEG_file;
    cfg.EEG_filepath = tEEG_filepath;
    clearvars tEEG_file tEEG_filepath
    headerICA.bad = bad;
    headerICA.badchans = badAll;
    headerICA.goodchans = setdiff(1:size(activations,1),headerICA.badchans);
    headerICA.goodchans = headerICA.goodchans(:)';
else
    badAll = [];
end

%--------------------------------------------------------------------------
% Write ICA activation (*_ICA.sef)
%--------------------------------------------------------------------------
disp('   Write results of ICA')
lab_write_sef(fullfile(cfg.ICA_filepath,cfg.ICA_file),activations,headerICA);

%--------------------------------------------------------------------------
% Write ICA parameter file (*_ICA.mat)
%--------------------------------------------------------------------------
cfg.ICA_param_file=[cfg.ICA_file(1:end-4) '.mat'];
save(fullfile(cfg.ICA_filepath,cfg.ICA_param_file),'W','cfg','header','ICAchans','headerICA');

%--------------------------------------------------------------------------
% Write exclude file (*_exclude_P.txt)
%--------------------------------------------------------------------------
Exclude_file=[cfg.ICA_file(1:end-4) '_exclude~.txt'];
fid=fopen(fullfile(cfg.ICA_filepath,Exclude_file),'w');
if ~isempty(badAll)
    fprintf(fid,num2str(badAll));
end
fclose(fid);

%--------------------------------------------------------------------------
% Write ICA activation topography (*_ICAtopo.sef)
%--------------------------------------------------------------------------
if isfield(headerICA,'events')
    headertopo = rmfield(headerICA,'events');
else
    headertopo = headerICA;
end
headertopo.numtimeframes = length(ICAchans);
cfg.ICAtopo_file=[cfg.ICA_file(1:end-4) 'topo.sef'];
lab_write_sef(fullfile(cfg.ICA_filepath,cfg.ICAtopo_file),topogood',headertopo);
ELS_file = fullfile(cfg.ICA_filepath,[cfg.ICAtopo_file(1:end-4) '.els']);
ELS_file = lab_write_locs(ELS_file,headertopo);
if ~isempty(ELS_file)
    % Write *.LM-file
    fidout=fopen(fullfile(cfg.ICA_filepath,[cfg.ICAtopo_file(1:end-4) '.lm']),'w');
    fprintf(fidout,[cfg.ICAtopo_file native2unicode([13 10])]);
    fprintf(fidout,[ELS_file native2unicode([13 10])]);
    fclose(fidout);
end

%--------------------------------------------------------------------------
% Write ICA activation topography maps (Topo*/*_ch.jpg)
%--------------------------------------------------------------------------
if ~isfield(cfg.ICA,'plottopo') | cfg.ICA.plottopo == true
    if isfield(headertopo,'locs')
        if max(headertopo.locs.z) > 0 | min(headertopo.locs.z) < 0
            settings.LOCS = headertopo.locs;
            settings.LOCS  = lab_locs2sph(settings.LOCS );
        end
    elseif isfield(cfg,'settings_path') & exist(fullfile(cfg.settings_path,'electrodes.els'),'file')
        settings.LOCS = lab_read_locs(fullfile(cfg.settings_path,'electrodes.els'));
    elseif isfield(cfg,'settings_path') & exist(fullfile(cfg.settings_path,'electrodes.sfp'),'file')
        settings.LOCS = lab_read_locs(fullfile(cfg.settings_path,'electrodes.sfp'));
    end
    if exist('settings','var')
        if size(topogood,2) == size(settings.LOCS.x,2)
            topo = topogood;
            plotchans = 1:size(topo,2);
        else
            plotchans = ICAchans;
        end
        if size(topo,2) == size(settings.LOCS.x,2)
            disp('   Write ICA activation topography maps')
            % Check for header.splitchans (e.g. MEG data)
            if isstruct(settings.LOCS) & isfield(header,'splitchans')
                nfigs = size(header.splitchans,1);
                topotmp = topo;
                clearvars topo plotchans
                for j = 1:size(header.splitchans,1)
                    plotchanstmp = plotchans(plotchans >= header.splitchans(j,1));
                    plotchanstmp = plotchanstmp(plotchanstmp <= header.splitchans(j,2));
                    topo{1,j} = topotmp(:,plotchanstmp);
                    locfile{1,j}.x = settings.LOCS.x(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.y = settings.LOCS.y(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.z = settings.LOCS.z(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.labels = settings.LOCS.labels(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.sph_radius = settings.LOCS.sph_radius(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.sph_theta = settings.LOCS.sph_theta(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.sph_phi = settings.LOCS.sph_phi(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.theta = settings.LOCS.theta(1,plotchanstmp); %#ok<AGROW>
                    locfile{1,j}.radius = settings.LOCS.radius(1,plotchanstmp); %#ok<AGROW>
                end
                clearvars topotmp plotchanstmp
            else
                nfigs = 1;
                topotmp = topo;
                plotchanstmp = plotchans;
                clearvars topo plotchans
                topo{1,1} = topotmp(:,plotchanstmp);
                locfile{1,1}.x = settings.LOCS.x(1,plotchanstmp);
                locfile{1,1}.y = settings.LOCS.y(1,plotchanstmp);
                locfile{1,1}.z = settings.LOCS.z(1,plotchanstmp);
                locfile{1,1}.labels = settings.LOCS.labels(1,plotchanstmp);
                locfile{1,1}.sph_radius = settings.LOCS.sph_radius(1,plotchanstmp);
                locfile{1,1}.sph_theta = settings.LOCS.sph_theta(1,plotchanstmp);
                locfile{1,1}.sph_phi = settings.LOCS.sph_phi(1,plotchanstmp);
                locfile{1,1}.theta = settings.LOCS.theta(1,plotchanstmp);
                locfile{1,1}.radius = settings.LOCS.radius(1,plotchanstmp);
                plotchans{1,1} = plotchanstmp; %#ok<NASGU>
                clearvars topotmp plotchanstmp
            end
            
            % create output folder
            cfg.ICAtopoimg_file = cfg.ICA_file(1:end-4);
            warning off %#ok<WNOFF>
            cfg.ICAtopoimg_filepath = fullfile(cfg.ICA_filepath,['Topo' cfg.ICAtopoimg_file]);
            mkdir (cfg.ICAtopoimg_filepath);
            warning on %#ok<WNON>
            
            % plot data
            settings.PLOT_file = [];
            settings.close = 0;
            settings.AddPlot = 0;
            for i = 1:size(topo{1,1},1)
                settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
                for j = 1:nfigs
                    subplot(1,nfigs,j)
                    settings.LOCS = locfile{1,j};
                    range = max(abs(topo{1,j}(i,:)));
                    PLOT.MinValue = -range;
                    PLOT.MaxValue = range;
                    PLOT.Color = lab_create_cmap('bluered');
                    settings = lab_plot_chans(topo{1,j}(i,:)',PLOT,settings);
                end
                lab_print_figure([fullfile(cfg.ICAtopoimg_filepath,cfg.ICAtopoimg_file) '-' num2str(i) '.jpg'],settings.handleF);
                close(settings.handleF);
            end
            clearvars topo;
            settings.handleF = figure('Visible','off','Color',[1 1 1]);
            stem(sum(abs(W),1));
            set(gca,'Xlim',[1 size(W,2)])
            title('ICA weights');
            lab_print_figure([fullfile(cfg.ICAtopoimg_filepath,cfg.ICAtopoimg_file) '-Weights.jpg'],settings.handleF);
            close(settings.handleF);
        end
    end
    clearvars headertopo settings PLOT
end

%--------------------------------------------------------------------------
% Write verbose file (*_ICA.vrb)
%--------------------------------------------------------------------------
Verbose_file=[cfg.ICA_file(1:end-4) '.vrb'];
fid=fopen(fullfile(cfg.ICA_filepath,Verbose_file),'w');
fprintf(fid,'ICA decomposition\n');
fprintf(fid,datestr(now,0));
fprintf(fid,'\n');
fprintf(fid,'Files:\n');
fprintf(fid,['EEG input file: ' cfg.EEG_file]);
fprintf(fid,'\n\n');
fprintf(fid,['Channels included:\n' num2str(ICAchans)]);
fprintf(fid,'\n\n');
fprintf(fid,['ICA components file: ' cfg.ICA_file]);
fprintf(fid,'\n');
fprintf(fid,['ICA topography file: ' cfg.ICAtopo_file]);
fprintf(fid,'\n');
fprintf(fid,['ICA parameter file: ' cfg.ICA_param_file]);
fprintf(fid,'\n\n');
if isfield(header,'badchans') & ~isempty(header.badchans)
    fprintf(fid,'Bad channels: ');
    fprintf(fid,num2str(header.badchans));
    fprintf(fid,'\n\n');
end
if ~isfield(bad,'error') | bad.error == 1
    fprintf(fid,'Error detecting bad activations!!\n');
end
fprintf(fid,'\n');
fclose(fid);

return