% Calculate Microstates using Kmeans-Clustering | T-AAHC Clustering
% (optional only timeframes at gfp-peaks are clustered)
% 
% [data,header,cfg] = lab_microstates(data,header,cfg)
%
% written by F. Hatz 2014

function [data,header,cfg] = lab_microstates(data,header,cfg)

global MicroTmp
if isfield(cfg,'firstfile') & cfg.firstfile == true
    MicroTmp = [];
end
    

if ~exist('data','var')
    [data,header,cfg] = lab_read_data;
end
if ~exist('header','var')
    header = lab_create_header(data);
end
if ~exist('cfg','var')
    cfg = [];
end

if isfield(cfg,'SKIP') & isfield(cfg.SKIP,'MICROST') & cfg.SKIP.MICROST == true;
    return
end

disp('Microstates')
if ~exist('cfg','var') | ~isfield(cfg,'MICROST') | ~isfield(cfg.MICROST,'eegsource');
    cfg.MICROST = [];
    [cfg,skipprocessing] = lab_set_microstates(cfg,header);
    if skipprocessing == 1
        return
    end
end

if ~isfield(cfg.MICROST,'folder') | isempty(cfg.MICROST.folder)
    cfg.MICROST.folder = 'Microstates';
end
if isfield(cfg,'Output_filepath') & exist(cfg.Output_filepath,'dir')
    Micro_filepath = fullfile(cfg.Output_filepath,cfg.MICROST.folder);
    Micro_file = cfg.Output_file;
elseif isfield(cfg,'EEG_filepath') & exist(cfg.EEG_filepath,'dir')
    Micro_filepath = fullfile(cfg.EEG_filepath,cfg.MICROST.folder);
    Micro_file = cfg.EEG_file;
elseif isfield(header,'EEG_filepath') & exist(header.EEG_filepath,'dir')
    Micro_filepath = fullfile(header.EEG_filepath,cfg.MICROST.folder);
    Micro_file = header.EEG_file;
else
    Micro_filepath = '';
    Micro_file = '';
end
if ~isempty(Micro_filepath)
    warning off %#ok<WNOFF>
    mkdir(Micro_filepath)
    cd(Micro_filepath);
    warning on %#ok<WNON>
    [~,~,~,Micro_fileS] = lab_filename(Micro_file);
end

% Calculate reference and interpolate bad
if strcmp(cfg.MICROST.eegsource,'montage')
    if cfg.MICROST.interpolate == 1
        disp('   Interpolate bad channels')
        [data,header] = lab_interpolate_bad(data,header);
    end
    if ~isfield(cfg.MICROST,'montage') | isempty(cfg.MICROST.montage)
        disp('   no montage loaded, use input structure')
        cfg.MICROST.montage = lab_create_montage(size(data,1),header);
    end
    if cfg.MICROST.montage(1,1).numchans ~= header.numchannels
        cfg.MICROST.montage = lab_reduce_montage(cfg.MICROST.montage,cfg,header,true);
    end
    if cfg.MICROST.montage(1,1).numchans >= header.numchannels
        disp('   invalid montage loaded, use input structure')
        montage = lab_create_montage(size(data,1),header);
    else
        montage = cfg.MICROST.montage;
    end
    [data,header] = lab_references(data,header,montage,cfg.MICROST);
elseif strcmp(cfg.MICROST.eegsource,'mean') | strcmp(cfg.MICROST.eegsource,'median') | strcmp(cfg.MICROST.eegsource,'laplacian')
    [data,header] = lab_references(data,header,cfg.MICROST.eegsource,cfg.MICROST);
    if cfg.MICROST.interpolate == 1
        disp('   Interpolate bad channels')
        [data,header] = lab_interpolate_bad(data,header);
    end
elseif strcmp(cfg.MICROST.eegsource,'chanref')
    [data,header] = lab_references(data,header,cfg.MICROST.eegsource,cfg.MICROST);
    if cfg.MICROST.interpolate == 1
        disp('   Interpolate bad channels')
        [data,header] = lab_interpolate_bad(data,header);
    end
elseif cfg.MICROST.interpolate == 1
    disp('   Interpolate bad channels')
    [data,header] = lab_interpolate_bad(data,header);
end

% reduce to data channels
if isfield(header,'numdatachannels')
    [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);
end

% exclude markers if necessary
if isfield (cfg.MICROST,'markerexclude') & iscell(cfg.MICROST.markerexclude) & ...
        ~isempty(cfg.MICROST.markerexclude) & isfield(header,'events') & ~isempty(header.events)
    disp(['   Exclude Markers: '  sprintf('%s ',cfg.MICROST.markerexclude{:})])
    [datatmp,headertmp,Idx] = lab_exclude_markers(data,header,cfg.MICROST.markerexclude);
    IdxBad = setdiff(1:size(data,2),Idx);
else
    datatmp = data;
    headertmp = header;
    Idx = 1:size(data,2);
    IdxBad = [];
end

if cfg.MICROST.dogfp == true
    disp('   Reduce to GFP-peaks for analysis')
    [~,~,gfp] = lab_compute_gfp(datatmp,headertmp);
    IdxHigh = find(diff(sign(diff(gfp)))== -2) + 1;
    IdxLow = (find(diff(sign(diff(gfp)))== 2) + 1);
    if min(IdxLow) > min(IdxHigh)
        IdxLow = [0 IdxLow];
    end
    if max(IdxLow) < max(IdxHigh)
        IdxLow = [IdxLow size(datatmp,2)];
    end
    
    [datatmp,headertmp] = lab_data2timeframes(datatmp,headertmp,IdxHigh);
    datatmp = datatmp / mean(gfp);
    IdxLow(2:end) = Idx(IdxLow(2:end));
    Idx = Idx(IdxHigh); %#ok<NASGU>
else
    IdxLow = [0 1:size(datatmp,2)];
end

if isempty(datatmp) | size(datatmp,2) <= 2
    disp('   Abort: no data to process')
    return
end

if isfield(cfg.MICROST,'doallfiles') & cfg.MICROST.doallfiles == true & isfield(cfg,'lastfile')
    % Collect files and process after collecting last file
    if ~isempty(MicroTmp) & isfield(MicroTmp,'data')
        if size(MicroTmp.data,1) == size(datatmp,1)
            [datatmp,headertmp] = lab_stitch(MicroTmp.data,MicroTmp.header,datatmp,headertmp,0);
            if isfield(headertmp,'locs') & isfield(MicroTmp.header,'locs')
                headertmp.locs.x = (headertmp.locs.x + MicroTmp.Nr*MicroTmp.header.locs.x) / (MicroTmp.Nr + 1);
                headertmp.locs.y = (headertmp.locs.y + MicroTmp.Nr*MicroTmp.header.locs.y) / (MicroTmp.Nr + 1);
                headertmp.locs.z = (headertmp.locs.z + MicroTmp.Nr*MicroTmp.header.locs.z) / (MicroTmp.Nr + 1);
                headertmp.locs = lab_locs2sph(headertmp.locs);
            end
            AllNr = MicroTmp.Nr + 1;
        else
            disp(['   Skip ' Micro_file '- not matching to previous files'])
            return
        end
    else
        AllNr = 1;
    end
    if cfg.lastfile == 0
        disp('   Keep data in memory and process next file')
        MicroTmp.data = datatmp;
        MicroTmp.header = headertmp;
        MicroTmp.Nr = AllNr;
        return
    else
        if isfield(cfg,'settings_path')
            Micro_filepath = fullfile(cfg.settings_path,cfg.MICROST.folder);
            warning off %#ok<WNOFF>
            mkdir(Micro_filepath)
            cd(Micro_filepath);
            warning on %#ok<WNON>
            Micro_file = ['Template_' num2str(AllNr) '-Files.sef'];
            [~,~,~,Micro_fileS] = lab_filename(Micro_file);
            lab_write_sef(fullfile(Micro_filepath,Micro_file),datatmp,headertmp);
        end
        disp(['   Calculate Microstates for ' Micro_file])
        data = datatmp;
        header = headertmp;
        skipdata = 1;
        MicroTmp = [];
    end
    clearvars AllNr
else
    skipdata = 0;
end

if isfield(cfg.MICROST,'template') & ~isempty(cfg.MICROST.template)
    % Do Fitting
    Template = cfg.MICROST.template;
    if size(Template,1) ~= size(datatmp,1)
        disp('   Abort: Different number of channels in Microstate-Template & Input-file')
        return
    end
    ClusterNr = size(Template,2);
    disp(['   Fit to Microstates-Template (' num2str(ClusterNr) ' templates)'])
    numclusters = 1;
    if isfield(cfg.MICROST,'AddZeros') & cfg.MICROST.AddZeros == true
        Template(:,end+1) = 0;
    end
    Result.CorrAll{1} = corr(Template,datatmp).^2;
    [CORR,Result.Clusters] = max(Result.CorrAll{1},[],1);
    Result.CORR = CORR;
    if isfield(cfg.MICROST,'AddZeros') & cfg.MICROST.AddZeros == true
        Result.Clusters(Result.Clusters == ClusterNr + 1) = 0;
        Template = Template(:,1:end-1);
    end
    Result.Nr = ClusterNr;
    Result.Template{1} = Template;
    clearvars Template
    
    flag_clustering = false;
else
    % Do Clustering
    if cfg.MICROST.minclusters < 2 | cfg.MICROST.minclusters > cfg.MICROST.maxclusters
        disp('   Abort: no valid maximal number of clusters')
        return
    end
    if isfield(cfg.MICROST,'method') & strcmp(cfg.MICROST.method,'T-AAHC')
        Result = lab_taahc(datatmp,cfg.MICROST);
    else
        Result = lab_kmeans(datatmp,cfg.MICROST);
    end
    ClusterNr = Result.Nr;
    numclusters = length(ClusterNr);

    flag_clustering = true;
end

% Convert GFP-data to full data
if (length(IdxLow) - 1) < size(data,2) & skipdata == 0
    Clusters = Result.Clusters;
    CORR = Result.CORR;
    if isfield(Result,'CorrAll')
        CorrAll = Result.CorrAll;
    end
    Result.Clusters = zeros(numclusters,size(data,2));
    Result.CORR = zeros(numclusters,size(data,2));
    for j = 1:numclusters;
        if isfield(Result,'CorrAll')
            Result.CorrAll{j} = zeros(size(CorrAll{j},1),size(data,2));
        end
        for i = 1:size(Clusters,2);
            Result.Clusters(j,IdxLow(i)+1:IdxLow(i+1)) = Clusters(j,i);
            Result.CORR(j,IdxLow(i)+1:IdxLow(i+1)) = CORR(j,i);
            if isfield(Result,'CorrAll')
                for m = 1:size(CorrAll{j},1)
                    Result.CorrAll{j}(m,IdxLow(i)+1:IdxLow(i+1)) = CorrAll{j}(m,i);
                end
            end
        end
    end
    clearvars Clusters i j m
    if ~isempty(IdxBad)
        Result.Clusters(:,IdxBad) = 0;
    end
end

% Correct for short segments
Result = reject_segments(Result,cfg.MICROST.rejectsegm);

% correct microstates orientation
for j = 1:length(Result.Template)
    Result.Template{j} = lab_orient_microstates(Result.Template{j},header);
end

% Calculate Stats
if size(Result.Clusters,2) == size(data,2)
    [Result,Stats] = lab_clustering_stats(Result,data);
else
    Stats = [];
end

% Create channels with Microstates-information
if size(Result.Clusters,2) == size(data,2) & skipdata == 0
    if ~isfield(header,'numdatachannels')
        header.numdatachannels = size(data,1);
    end
    channels = cellstr(header.channels);
    data = cat(1,data,Result.Clusters);
    for i = 1:numclusters
        channels{end+1,1} = ['Micro' num2str(ClusterNr(i))]; %#ok<AGROW>
    end
    for i = 1:numclusters
        data = cat(1,data,Stats(i).CORR);
        channels{end+1,1} = ['Corr' num2str(ClusterNr(i))]; %#ok<AGROW>
    end
    header.channels = char(channels);
    header.numchannels = size(data,1);
    header.numauxchannels = header.numchannels - header.numdatachannels;
    
    % Create events for microstates
    if length(ClusterNr) > 1
        if isfield(cfg.MICROST,'doevents_mode') & strcmp(cfg.MICROST.doevents_mode,'auto')
            warning off %#ok<WNOFF>
            tmp = corner(log(Result.KL),1:length(Result.KL));
            if ~isempty(tmp);
                header.microstates = ClusterNr(tmp(1));
            else
                header.microstates = [];
            end
            warning on %#ok<WNON>
        elseif isfield(cfg.MICROST,'doevents_mode') & strcmp(cfg.MICROST.doevents_mode,'auto first')
            if isempty(cfg.MICROST.doevents) | cfg.MICROST.doevents < 2
                warning off %#ok<WNOFF>
                tmp = corner(log(Result.KL),1:length(Result.KL));
                if ~isempty(tmp);
                    cfg.MICROST.doevents = ClusterNr(tmp(1));
                else
                    cfg.MICROST.doevents = [];
                end
                warning on %#ok<WNON>
            end
            header.microstates = cfg.MICROST.doevents;
        elseif isfield(cfg.MICROST,'doevents') & cfg.MICROST.doevents > 0 & cfg.MICROST.doevents <= max(ClusterNr)
            header.microstates = cfg.MICROST.doevents;
        else
            header.microstates = [];
        end
    elseif length(ClusterNr) == 1
        header.microstates = ClusterNr;
    end
    if ~isempty(header.microstates)
        [data,header] = lab_microstates2events(data,header,header.microstates);
    end
end

% Write Results
if isempty(Micro_filepath)
    return
end

% Write Mat-Container
save(fullfile(Micro_filepath,[Micro_fileS '_Microstates.mat']),'Result','Stats','header','-v7.3');

% Write Microstates Topos
if flag_clustering == true
    Topos = Result.Template;
    Label = 'MicrostatesTopo';
else
    Topos{1} = Stats(1).ClusterTopo;
    Label = 'Clusters';
end

settings = [];
for i = 1:numclusters
    % create output folder for topos
    Output_filepath2 = fullfile(Micro_filepath,[num2str(ClusterNr(i)) Label]);
    warning off %#ok<WNOFF>
    mkdir(Output_filepath2)
    warning on %#ok<WNON>
    
    % write topo as Cartool.sef-file
    headertopo = header;
    headertopo.numtimeframes = ClusterNr(i);
    lab_write_sef(fullfile(Output_filepath2,[Micro_fileS '_' num2str(ClusterNr(i)) 'Clusters.sef']),Topos{i},headertopo);
    
    if isfield(headertopo,'locs') & size(Topos{i},1) == size(headertopo.locs.x,2) & ...
            (max(headertopo.locs.z) > 0 | min(headertopo.locs.z) < 0)
        if strcmp(cfg.MICROST.doplots,'2D')
            disp(['   Write Microstate Topos (' num2str(ClusterNr(i)) ' Clusters)'])
            if ~isfield(settings,'LOCS')
                settings.LOCS = headertopo.locs;
                settings.LOCS = lab_locs2sph(settings.LOCS);
                settings.PLOT_file = '';
            end
            PLOT.close = 0;
            PLOT.AddPlot = 0;
            for j = 1:size(Topos{i},2)
                settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
                range = max(abs(Topos{i}(:,j)));
                PLOT.MinValue = -range;
                PLOT.MaxValue = range;
                PLOT.Color = lab_create_cmap('bluered');
                settings = lab_plot_chans(Topos{i}(:,j),PLOT,settings);
                lab_print_figure(fullfile(Output_filepath2,['Cluster_' num2str(j) '.jpg']),settings.handleF);
                close(settings.handleF);
            end
            clearvars PLOT range f
        elseif strcmp(cfg.MICROST.doplots,'3D')
            Cmap = lab_create_cmap('bluered');
            PLOT.facecolor = [];
            PLOT.alpha = 1;
            PLOT.sizedots = 4;
            PLOT.plotlabels = true;
            PLOT.labelsize = 16;
            PLOT.labeldistance = 1;
            PLOT.plotdots = true;
            PLOT.invisible = true;
            for j = 1:size(Topos{i},2)
                Vars = (Topos{i}(:,j)- min(Topos{i}(:,j))) / (max(Topos{i}(:,j)) - min(Topos{i}(:,j)));
                for V = 1:length(Vars)
                    PLOT.dotscolor(V,:) = Cmap(round(Vars(V)*63+1),:);
                end
                Fig = lab_plot_mesh(headertopo.locs,PLOT);
                view(gca,[0 90]);
                lab_print_figure(fullfile(Output_filepath2,['Cluster_' num2str(j) '.jpg']),Fig);
                close(Fig);
            end
            clearvars Cmap Vars PLOT Fig
        end
    end
end

if isempty(Stats)
    return
end

for i = 1:numclusters
    Output_filepath2 = fullfile(Micro_filepath,[num2str(ClusterNr(i)) Label]);
    warning off %#ok<WNOFF>
    mkdir(Output_filepath2)
    warning on %#ok<WNON>
    
    xlsout = {'Cluster';'Mean EV';'Mean Correlation';'Max Correlation';'Max Correlation Index'; ...
        'Max GFP';'Max GFP Index';'Mean Duration';'Minimal Duration'; ...
        'Maximal Duration';'Std Duration';'Sum Duration'};
    xlsout = cat(2,xlsout,cat(1,num2cell(1:Stats(i).Nr),num2cell(Stats(i).MeanEV(:)'), ...
        num2cell(Stats(i).MeanCorr(:)'),num2cell(Stats(i).MaxCorr(:)'),num2cell(Stats(i).MaxCorrIdx(:)'), ...
        num2cell(Stats(i).MaxGfp(:)'),num2cell(Stats(i).MaxGfpIdx(:)'), ...
        num2cell(Stats(i).MeanDuration(:)'),num2cell(Stats(i).MinDuration(:)'), ...
        num2cell(Stats(i).MaxDuration(:)'),num2cell(Stats(i).StdDuration(:)'), ...
        num2cell(Stats(i).SumDuration(:)')));
    xlstmp = {};
    for j = 1:Stats(i).Nr
        xlstmp{j,1} = ['NextCluster' num2str(j)]; %#ok<AGROW>
    end
    xlstmp = cat(2,xlstmp,num2cell(Stats(i).NextCluster ./ repmat(sum(Stats(i).NextCluster,1),Stats(i).Nr,1)));
    xlsout = cat(1,xlsout,xlstmp);
    xlstmp = {};
    for j = 1:Stats(i).Nr
        xlstmp{j,1} = ['PreviousCluster' num2str(j)]; %#ok<AGROW>
    end
    xlstmp = cat(2,xlstmp,num2cell(Stats(i).PreviousCluster ./ repmat(sum(Stats(i).PreviousCluster,1),Stats(i).Nr,1)));
    xlsout = cat(1,xlsout,xlstmp);
    lab_write_xls(fullfile(Output_filepath2,[Micro_fileS '_' num2str(ClusterNr(i)) 'Clusters.xlsx']),xlsout);
    clearvars xlsout xlstmp j
end

% Write GEV and Variance of clustering
Fnames = fieldnames(Stats);
Idx = [find(strcmp(Fnames,'Nr')) find(strcmp(Fnames,'GEV')) find(strcmp(Fnames,'MeanVariance')) ...
    find(strcmp(Fnames,'KL')) find(strcmp(Fnames,'CV'))];
Stats2 = struct2cell(Stats);
Stats2 = permute(Stats2(Idx,1,:),[1 3 2]);
xlsout = cat(1,{'Number of clusters','GEV','Mean Variance','Krzanowski-Lai','Cross-Validation'},Stats2');
if isfield(Result,'Niter')
    xlsout = cat(2,xlsout,cat(1,{'Iterations'},num2cell(Result.Niter(:))));
end
lab_write_xls(fullfile(Micro_filepath,'Clustering.xlsx'),xlsout);

% Plot GEV and Variance of clustering
if numclusters > 1
    Stats2 = cell2mat(Stats2);
    ClusterNr = Stats2(1,:);
    fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Clustering -- GEV','Visible','off');
    
    plot(Stats2(2,:));
    title('Clustering -- Global Explained Variance')
    set(gca,'XTick',1:length(ClusterNr),'XTickLabel',cellstr(num2str(ClusterNr')));
    xlabel('Number of clusters');
    lab_print_figure(fullfile(Micro_filepath,'Microstates-GEV.tif'),fig1);
    
    plot(Stats2(3,:));
    title('Clustering -- Mean Variance')
    set(gca,'XTick',1:length(ClusterNr),'XTickLabel',cellstr(num2str(ClusterNr')));
    xlabel('Number of clusters');
    lab_print_figure(fullfile(Micro_filepath,'Microstates-MeanVar.tif'),fig1);
    
    plot(Stats2(4,:));
    title('Clustering -- Krzanowski-Lai')
    set(gca,'XTick',1:length(ClusterNr),'XTickLabel',cellstr(num2str(ClusterNr')));
    xlabel('Number of clusters');
    lab_print_figure(fullfile(Micro_filepath,'Microstates-KrzanowskiLai.tif'),fig1);
    
    plot(Stats2(5,:));
    title('Clustering -- Cross Validation')
    set(gca,'XTick',1:length(ClusterNr),'XTickLabel',cellstr(num2str(ClusterNr')));
    xlabel('Number of clusters');
    lab_print_figure(fullfile(Micro_filepath,'Microstates-CrossValidation.tif'),fig1);
    
    close(fig1);
end

% Write EEG/MEG with microstates
if isfield(header,'microstates') & ~isempty(header.microstates)
    lab_write_sef(fullfile(Micro_filepath,[Micro_fileS '_' num2str(header.microstates) 'Microstates.sef']),data,header);
end

% do calculations on Microstates
Outputfile = cfg.Output_file;
OutputfileS = cfg.Output_fileS;
Output_filepath = cfg.Output_filepath;
cfg.Output_file = Micro_file;
cfg.Output_fileS = Micro_fileS;
cfg.Output_filepath = Micro_filepath;
cfg.SKIP.MICROST = true;
cfg = lab_postprocessing(data,header,cfg,true);
cfg.Output_file = Outputfile;
cfg.Output_fileS = OutputfileS;
cfg.Output_filepath = Output_filepath;
clearvars Outputfile OutputfileS

end

function Result = reject_segments(Result,rejectsegm)
    if rejectsegm > 1
        fprintf(['   Correct for short segments (<' num2str(rejectsegm) ' TF)'])
        % for i = 1:size(Result.Clusters,1)
        for i = 1:size(Result.Clusters,1)
            Clustertmp = Result.Clusters(i,:);
            Clustertmp(Clustertmp==0) = Result.Nr(i) + 1;
            CORR = Result.CORR(i,:);
            CorrAll = Result.CorrAll{i};
            Idx = find(abs(diff(Clustertmp))>0);
            Idx = [0 Idx size(Clustertmp,2)]; %#ok<AGROW>
            L = Idx(2:end) - Idx(1:end-1);
            IdxL = find(L <= rejectsegm);
            while ~isempty(IdxL)
                for j = 1:length(IdxL)
                    T = [];
                    if Idx(IdxL(j)) > 0
                        T = union(T,Clustertmp(Idx(IdxL(j))));
                    end
                    if Idx(IdxL(j)+1)+1 <= length(Clustertmp)
                        T = union(T,Clustertmp(Idx(IdxL(j)+1)+1));
                    end
                    T = setdiff(T,Result.Nr(i) + 1);
                    tmp = mean(CorrAll(T,Idx(IdxL(j))+1:Idx(IdxL(j)+1)),2);
                    [~,tmp] = max(tmp);
                    T = T(tmp);
                    Clustertmp(Idx(IdxL(j))+1:Idx(IdxL(j)+1)) = T;
                    CORR(Idx(IdxL(j))+1:Idx(IdxL(j)+1)) = CorrAll(T,Idx(IdxL(j))+1:Idx(IdxL(j)+1));
                end
                Idx = find(abs(diff(Clustertmp))>0);
                Idx = [0 Idx size(Clustertmp,2)]; %#ok<AGROW>
                L = Idx(2:end) - Idx(1:end-1);
                IdxL = find(L <= rejectsegm);
            end
            Clustertmp(Clustertmp==Result.Nr(i) + 1) = 0;
            Result.Clusters(i,:) = Clustertmp;
            Result.CORR(i,:) = CORR;
            fprintf('.');
        end
        disp(':');
    end
end