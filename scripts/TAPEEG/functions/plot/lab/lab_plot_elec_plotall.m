function lab_plot_elec_plotall(cfg)

Nchans = size(cfg.LOCS.x,2);
if isfield(cfg,'Mappings') & ~isempty(cfg.Mappings)
    Nchans = size(cfg.Mappings.mappings,2);
end
[DATA,PLOT,PlotList] = lab_plot_prepare_data(cfg,[1 1 1]);
if isempty(DATA)
    disp(['Wrong input data, must be exactly ' num2str(Nchans) ' values'])
    return
end
[DATA2,PLOT2] = lab_plot_mix_data(DATA,PLOT,PlotList);
[DATA2,PLOT2] = lab_plot_create_title(DATA2,PLOT2);

for plotnr = 1:size(DATA2,1)
    Modes = fieldnames(DATA2(plotnr));
    FlagM = false;
    Ptitle = '';
    settings = cfg;
    Legend = {};
    for i = 1:length(Modes)
        if isfield(DATA2(plotnr).(Modes{i}),'data') & ~isempty(DATA2(plotnr).(Modes{i}).data) & ...
                (~isfield(PLOT2(plotnr).(Modes{i}),'Valid') | PLOT2(plotnr).(Modes{i}).Valid == true)
            if FlagM == false
                PLOT2(plotnr).(Modes{i}).AddPlot = false;
                FlagM = true;
            else
                PLOT2(plotnr).(Modes{i}).AddPlot = true;
            end
            if strcmp(Modes{i},'Nodes')
                data = DATA2(plotnr).(Modes{i}).data;
                plotslider = false;
                plotsingle = true;
            elseif strcmp(Modes{i},'Connections')
                data = lab_tril2matrix(DATA2(plotnr).(Modes{i}).data,Nchans);
                if size(data,3) > 1
                    data = mean(data,3);
                    disp('More than one matrix for single plot selected, show average matrix')
                end
                data(1:size(data,1)+1:end) = NaN;
                plotslider = true;
                plotsingle = true;
            elseif strcmp(Modes{i},'Surface')
                data = DATA2(plotnr).(Modes{i}).data;
                plotslider = false;
                plotsingle = false;
            else
                data = [];
            end
            if ~isempty(data) & length(data) == Nchans
                if isfield(PLOT2(plotnr),'Name')
                    Ptitle =  PLOT2(plotnr).Name;
                else
                    Ptitle = 'Plot';
                end
                PLOT2(plotnr).(Modes{i}).Title = Ptitle;
                if ~isfield(settings,'handleF') | isempty(settings.handleF) | ~ishandle(settings.handleF)
                    PLOT2(plotnr).(Modes{i}).AddPlot = false;
                end
                if isfield(PLOT2(plotnr).(Modes{i}),'Legend') & ~isempty(PLOT2(plotnr).(Modes{i}).Legend)
                    Legend{end+1} = PLOT2(plotnr).(Modes{i}).Legend; %#ok<AGROW>
                end
                if ~isfield(cfg,'onlyhead') | cfg.onlyhead == false
                    if cfg.plothead == true
                        plothead = 1;
                    else
                        plothead = 0;
                    end
                else
                    plothead = 2;
                end
                settings = lab_plot_chans(data,PLOT2(plotnr).(Modes{i}),settings,plotsingle,plotslider,plothead,Legend);
            end
        end
    end
    if isfield(settings,'handleF') & ~isempty(settings.handleF) & ishandle(settings.handleF)
        if cfg.Store == true
            [~,DATA_filepath,DATA_format,DATA_fileS] = lab_filename(cfg.DATA_file);
            DATA_fileS = [DATA_fileS '_' Ptitle]; %#ok<AGROW>
            lab_print_figure(fullfile(DATA_filepath,[DATA_fileS '.' DATA_format]),settings.handleF);
        end
    else
        disp('No valid data for printing')
    end
end

end