% Helper file for lab_plot_IS
%
% Written by F. Hatz 2013

function [nodes,surface,cfg] = lab_plot_create_nodes(matrix,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('matrix','var') | isempty(matrix)
    nodes = [];
    surface = [];
    return
end

if isstruct(matrix) & isfield(matrix,'matrix')
    matrix = matrix.matrix;
elseif ~isnumeric(matrix) | size(matrix,1) ~= size(matrix,2)
    nodes = [];
    surface = [];
    return
end
nodes = matrix(1:(size(matrix,1)+1):end);
surface = zeros(1,size(nodes,2));
matrix(1:size(matrix,1)+1:end) = 0;

if ~isfield(cfg,'nodemethod')
    cfg.nodemethod = 'No nodes';
    if ~isfield(cfg,'nodenormalize')
        cfg.nodenormalize = false;
    end
    cfg.surfacemethod = 'No surface';
    if ~isfield(cfg,'surfacenormalize')
        cfg.surfacenormalize = false;
    end
    
    Prompt = cell(0,2);
    Formats = {};
    
    Prompt(end+1,:) = {'Nodes','nodemethod'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    if max(nodes) == 0 | min(nodes) == 1
        Formats(end,1).items = {'No nodes','Degree','Betweenness centrality','Eigenvector centrality', ...
            'Clustering coefficient','File','Small dots','Large dots'};
    else
        Formats(end,1).items = {'No nodes','Degree','Betweenness centrality','Eigenvector centrality', ...
            'Clustering coefficient','File','Small dots','Large dots','Diagonal data'};
    end
    
    Prompt(end+1,:) = {'Normalize','nodenormalize'};
    Formats(end+1,1).type = 'check';
    
    if isfield(cfg,'dosurface') & cfg.dosurface == true
        Prompt(end+1,:) = {'Surface','surfacemethod'};
        Formats(end+1,1).type = 'list';
        Formats(end,1).style = 'popupmenu';
        Formats(end,1).format = 'input';
        Formats(end,1).items = {'No surface','Betweenness centrality','Eigenvector centrality', ...
            'Clustering coefficient','File'};
        
        Prompt(end+1,:) = {'Normalize','surfacenormalize'};
        Formats(end+1,1).type = 'check';
    end
    
    [cfg,Cancelled] = inputsdlg(Prompt,'Nodes',Formats,cfg,2);
    if Cancelled == 1
        nodes = surface;
        return
    end
    pause(0.2);
end

switch cfg.nodemethod
    case 'Diagonal data'
        disp('    take diagonal data of matrix for nodes information')
    case 'Degree'
        if size(unique(matrix),1) == 2
            nodes = degrees_und(matrix);
        else
            nodes = mean(matrix,1);
        end
    case 'Betweenness centrality'
        disp('    calculate betweenness centrality')
        if size(unique(matrix),1) == 2
            [~,nodes] = edge_betweenness_bin(matrix);
            nodes = nodes' / ((size(matrix,1)-1)*(size(matrix,1)-2));
        else
            matrixI = -abs(matrix);
            matrixI(matrixI == 0) = NaN;
            matrixI(1:(size(matrix,1)+1):end) = NaN;
            matrixI = (matrixI - min(matrixI(:)));
            matrixI = matrixI / max(matrixI(:));
            matrixI(1:(size(matrix,1)+1):end) = 0;
            matrixI(isnan(matrixI)) = Inf;
            [~,nodes]=edge_betweenness_wei(matrixI);
            nodes = nodes' ./ ((size(matrix,1)-1)*(size(matrix,1)-2));
            clearvars matrixI
        end
    case 'Eigenvector centrality'
        disp('    calculate eigenvector centrality')
        nodes = eigenvector_centrality_und(matrix)';
    case 'Clustering coefficient'
        disp('    calculate Clustering coefficients')
        if size(unique(matrix),1) == 2
            nodes = clustering_coef_bu(abs(matrix))';
        else
            nodes = clustering_coef_wu(abs(matrix))';
        end
    case 'File'
        for i = 1:size(nodes,2)
            labels{i,1} = ['N_' num2str(i,'%03g')]; %#ok<AGROW>
        end
        clearvars i
        [nodes,~,cfgnode,skipprocessing] = lab_read_plotdata([],labels,cfg,1);
        if skipprocessing == 0
            varnames = cfgnode.varnames;
            clearvars cfgnode
            if size(nodes,2) == 1
                nodes = nodes';
            end
            if size(nodes,1) > 1
                selnode = listdlg('PromptString','Variable for nodes','SelectionMode', ...
                    'single','ListString',varnames);
                if ~isempty(selnode)
                    nodes = nodes(selnode,:);
                    tmp = setdiff(1:size(varnames,2),selnode);
                    Svarnames = varnames(tmp);
                    Surface = nodes(tmp,:);
                    clearvars tmp selnode varnames
                else
                    nodes = zeros(1,size(matrix,1));
                end
            end
        end
        if size(nodes,2) ~= size(labels,1)
            nodes = zeros(1,size(matrix,1));
            disp('  number of nodes not matching, replaced by zeros')
        end
    case 'Small dots'
        nodes = ones(1,size(matrix,1))*0.5;
    case 'Large dots'
        nodes = ones(1,size(matrix,1));
    otherwise
        nodes = zeros(1,size(matrix,1));
end
if cfg.nodenormalize == true
    mx = max(nodes);
    mn = min(nodes);
    if mx > mn & mn > 0
        nodes = (nodes / mx);
    elseif mx > mn
        nodes = ((nodes-mn) / (mx-mn));
    elseif mx ~= mn
        nodes = zeros(1,size(matrix,1));
    end
    clearvars mx mn
end

switch cfg.surfacemethod
    case 'Betweenness centrality'
        disp('    calculate betweenness centrality')
        if size(unique(matrix),1) == 2
            [~,surface] = edge_betweenness_bin(matrix);
            surface = surface' / ((size(matrix,1)-1)*(size(matrix,1)-2));
        else
            matrixI = -abs(matrix);
            matrixI(matrixI == 0) = NaN;
            matrixI(1:(size(matrix,1)+1):end) = NaN;
            matrixI = (matrixI - min(matrixI(:)));
            matrixI = matrixI / max(matrixI(:));
            matrixI(1:(size(matrix,1)+1):end) = 0;
            matrixI(isnan(matrixI)) = Inf;
            [~,surface]=edge_betweenness_wei(matrixI);
            surface = surface' ./ ((size(matrix,1)-1)*(size(matrix,1)-2));
            clearvars matrixI
        end
    case 'Eigenvector centrality'
        disp('    calculate eigenvector centrality')
        surface = eigenvector_centrality_und(matrix)';
    case 'Clustering coefficient'
        disp('    calculate Clustering coefficients')
        if size(unique(matrix),1) == 2
            surface = clustering_coef_bu(abs(matrix))';
        else
            surface = clustering_coef_wu(abs(matrix))';
        end
    case 'File'
        if exist('Surface','var') & exist('Svarnames','var')
            surface = Surface;
            varnames = Svarnames;
        else
            for i = 1:size(surface,2)
                labels{i,1} = ['N_' num2str(i,'%03g')]; %#ok<AGROW>
            end
            clearvars i
            [surface,~,cfgsurface,skipprocessing] = lab_read_plotdata([],labels,cfg,1);
            if skipprocessing == 1
                surface = zeros(1,size(matrix,1));
                return
            end
            varnames = cfgsurface.varnames;
            clearvars cfgnode
        end
        if size(surface,2) == 1
            surface = surface';
        end
        if size(surface,1) > 1 & size(surface,2) == size(labels,1)
            selsurf = listdlg('PromptString','Variable for surface','SelectionMode', ...
                'single','ListString',varnames);
            if ~isempty(selsurf)
                surface = surface(selsurf,:);
                varnames = varnames(selsurf);
            else
                surface = zeros(1,size(matrix,1));
                return
            end
        end
        PLOT.Name = regexprep(varnames{1},'_',' ');
        FORMAT{strcmp(fieldnames(PLOT),'Name')} = 'text';
        PLOT.MinValue = min(surface);
        PLOT.MaxValue = max(surface);
        if length(unique(surface)) == 2 & min(surface) == 0 & max(surface) == 1
            PLOT.Color = [0 1 0];
        else
            PLOT.Color = 'single color';
            FORMAT{strcmp(fieldnames(PLOT),'Color')} = {'single color','bluered','autumn','bone','colorcube','cool','copper','gray','hot','hsv','jet'};
        end
        PLOT.alpha = 0.2;
        PLOT.VolumePlot = false;
        PLOT = inputsdlg(PLOT,'Plot settings',FORMAT);
        if ischar(PLOT.Color) & strcmp(PLOT.Color,'single color')
            PLOT2.Color = [0 1 0];
            PLOT2 = inputsdlg(PLOT2,'Color',FORMAT);
            PLOT.Color = PLOT2.Color;
        end
        if abs(PLOT.MaxValue-PLOT.MinValue) > 0
            surface = ((surface-PLOT.MinValue) / (PLOT.MaxValue-PLOT.MinValue));
        else
            surface = zeros(1,size(matrix,1));
        end
        cfg.PLOT.MinValue = PLOT.MinValue;
        cfg.PLOT.MaxValue = PLOT.MaxValue;
        cfg.PLOT.alpha = [0.2 PLOT.alpha];
        cfg.PLOT.Color1 = PLOT.Color;
        cfg.PLOT.VolumePlot = PLOT.VolumePlot;
        cfg.PLOT.PlotAdd = false;
    otherwise
        surface = zeros(1,size(matrix,1));
end
mx = max(surface);
mn = min(surface);
if mx > mn & mn > 0 & cfg.surfacenormalize == false
    surface = (surface / mx);
elseif mx > mn
    surface = ((surface-mn) / (mx-mn));
elseif mx ~= mn
    surface = zeros(1,size(matrix,1));
end
clearvars mx mn

return