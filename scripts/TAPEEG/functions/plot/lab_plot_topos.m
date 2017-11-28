function lab_plot_topos

set.SEARCH.searchstring = {'*topo.sef'};
Files = lab_search_files(set);
Filelist = Files.Filelist;
if isempty(Filelist)
    return
end

for filenr = 1:length(Filelist)
    [data,header,cfg] = lab_read_data(Filelist{filenr});
    if ~isfield(header,'locs') | isempty(header.locs)
        disp('Abort, no LOCS information')
        return
    end
    
    settings.LOCS = header.locs;
    if isfield(header,'numdatachannels')
        [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);
    end
    
    if size(data,1) ~= size(settings.LOCS.x,2)
        disp('Abort, wrong LOCS information')
        return
    end
    
    disp('   Write topography maps')
    % Check for header.splitchans (e.g. MEG data)
    if isstruct(settings.LOCS) & isfield(header,'splitchans')
        nfigs = size(header.splitchans,1);
        clearvars topo plotchans
        for j = 1:size(header.splitchans,1)
            plotchanstmp = header.splitchans(j,1):size(data,1);
            plotchanstmp = plotchanstmp(plotchanstmp <= header.splitchans(j,2));
            topo{1,j} = data(plotchanstmp,:); %#ok<AGROW>
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
        clearvars topo plotchans
        topo{1,1} = data;
        locfile{1,1} = settings.LOCS;
        clearvars topotmp plotchanstmp
    end
    
    % create output folder
    cfg.ICAtopoimg_file = cfg.EEG_file(1:end-4);
    warning off %#ok<WNOFF>
    cfg.ICAtopoimg_filepath = fullfile(cfg.EEG_filepath,['Topo' cfg.ICAtopoimg_file]);
    mkdir (cfg.ICAtopoimg_filepath);
    warning on %#ok<WNON>
    
    % plot data
    settings.PLOT_file = [];
    settings.close = 0;
    settings.AddPlot = 0;
    for i = 1:size(topo{1,1},2)
        settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
        for j = 1:nfigs
            subplot(1,nfigs,j)
            settings.LOCS = locfile{1,j};
            range = max(abs(topo{1,j}(:,i)));
            PLOT.MinValue = -range;
            PLOT.MaxValue = range;
            PLOT.Color = lab_create_cmap('bluered');
            settings = lab_plot_chans(topo{1,j}(:,i),PLOT,settings);
        end
        lab_print_figure([fullfile(cfg.ICAtopoimg_filepath,cfg.ICAtopoimg_file) '-' num2str(i) '.jpg'],settings.handleF);
        close(settings.handleF);
    end
    clearvars topo;
    
    if isfield(header,'W')
        settings.handleF = figure('Visible','off','Color',[1 1 1]);
        stem(sum(abs(header.W),1));
        set(gca,'Xlim',[1 size(header.W,2)])
        title('ICA weights');
        lab_print_figure([fullfile(cfg.ICAtopoimg_filepath,cfg.ICAtopoimg_file) '-Weights.jpg'],settings.handleF);
        close(settings.handleF);
    end
end

clearvars header settings PLOT