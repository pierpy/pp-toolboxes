% Helper script for lab_calculate_connectivity
%
% written by F.Hatz 2014

function cfg = lab_calculate_connect_plotmatrix(Result,cfg)

disp('    Write matrices')

if isempty(Result)
    return
end
if ~isfield(cfg,'Output_file') | ~isfield(cfg,'Output_filepath')
    return
end

Connect_filepath = cfg.Output_filepath;
if isfield(cfg,'patient') & ~isempty(cfg.patient)
    Connect_file = [cfg.patient '_Connectivity'];
else
    Connect_file = 'Connectivity';
end

frequencies = fieldnames(Result);
tmp = [];
for i = 1:size(frequencies,1)
    if isstruct(Result.(frequencies{i,1})) & strcmp(frequencies{i,1}(1),'F')
        tmp = [tmp i];
    end
end
if isempty(tmp)
    frequencies = 0;
else
    frequencies = frequencies(tmp,1);
end

for Nfreq = 1:size(frequencies,1)
    if ~isnumeric(frequencies)
        result = Result.(frequencies{Nfreq,1});
        Outputname = frequencies{Nfreq,1};
    else
        result = Result;
        Outputname = '';
    end
    
    variables = fieldnames(result);
    tmp = [];
    for i = 1:size(variables,1)
        if isnumeric(result.(variables{i,1})) & size(result.(variables{i,1}),1) > 1 & ...
                size(result.(variables{i,1}),1) == size(result.(variables{i,1}),2) & ...
                ~strcmp(variables{i,1},'SL_hit') & ~strcmp(variables{i,1},'SLc_hit') & ...
                ~strcmp(variables{i,1},'wplv_angle')
            tmp = [tmp i];
        end
    end
    if isempty(tmp)
        return
    else
        variables = variables(tmp,1);
    end
    
    for Nvar = 1:size(variables,1)
        if size(result.(variables{Nvar,1}),3) <= 100
            if isfield(result,'timestamp') & ~isempty(result.timestamp) & length(result.timestamp) == size(result.(variables{Nvar,1}),3)
                timestamp = result.timestamp;
            else
                for i = 1:size(result.(variables{Nvar,1}),3)
                    timestamp{i,1} = num2str(i);
                end
            end
            for nmatrix = 1:size(result.(variables{Nvar,1}),3)
                matrix = result.(variables{Nvar,1})(:,:,nmatrix);
                matrix(1:size(matrix,1)+1:end) = 0;
                matrix_file = fullfile(Connect_filepath,[Connect_file '_' Outputname '_' timestamp{nmatrix} '_' upper(variables{Nvar,1}) '_matrix.tif']);
                if ~exist(matrix_file,'file')
                    cfg.CONNECT.PLOT = plot_matrix(matrix,matrix_file,cfg.CONNECT.PLOT);
                end
            end
        end
    end
end

end

function settings = plot_matrix(matrix,matrix_file,settings)
    disp(['    Plot matrix: ' matrix_file])
    settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
    settings.PLOT_file = [];
    settings.close = 0;
    if isfield(settings,'threshold') & ~isempty(settings.threshold)
        tmp = lab_rm_diagonal(matrix);
        tmp = sort(tmp(:));
        tmp2 = round((settings.threshold/100)*length(tmp));
        if tmp2 < 1; tmp2 = 1; end
        if tmp2 > length(tmp); tmp2 = length(tmp); end
        matrix(matrix < tmp(tmp2)) = 0;
    end
    matrix(1:size(matrix,1)+1:end) = 0;
    if max(matrix(:)) == 0 & min(matrix(:)) == 0
        PLOT.MinValueE = 0;
        PLOT.MaxValueE = 1;
    else
        PLOT.MinValueE = min(matrix(:));
        PLOT.MaxValueE = max(matrix(:));
    end
    PLOT.ColorE = lab_create_cmap(settings.ColorE);
    PLOT.SizeE = settings.SizeE;
    PLOT.Color = lab_create_cmap(settings.Color);
    PLOT.Size = settings.Size;
    if isfield(settings,'nodemethod')
        nodes = lab_plot_create_nodes(matrix,settings);
        matrix(1:size(matrix,1)+1:end) = nodes;
        if max(nodes(:)) == 0 & min(nodes(:)) == 0
            PLOT.MinValue = 0;
            PLOT.MaxValue = 1;
        else
            PLOT.MinValue = min(nodes(:));
            PLOT.MaxValue = max(nodes(:));
        end
    end
    settings = lab_plot_chans(matrix,PLOT,settings);
    lab_print_figure(matrix_file,settings.handleF);
    close(settings.handleF);
end