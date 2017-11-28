% Plot matrix
%
% lab_plot_matrix(matrix)
%
%       matrix     array (nchans x nchans)
%
% Written by F. Hatz 2013

function fig1 = lab_plot_matrix(matrix,doinvisible,nomenu)

fig1 = [];
header = [];
Graph = [];

if ~exist('nomenu','var')
    nomenu = false;
end
if ~exist('doinvisible','var')
    doinvisible = false;
end

if exist('matrix','var') & ischar(matrix)
    if ~exist(matrix,'file')
        clearvars matrix
    else
        [~,matrix_filepath,~,matrix_fileS] = lab_filename(matrix);
        matrix = lab_read_data(matrix);
        cd(matrix_filepath);
    end
end

if ~exist('matrix','var') | ~isnumeric(matrix)
    [matrix_file,matrix_filepath] = uigetfile('*.*','Select matrix-file','MultiSelect','on');
    if ~iscell(matrix_file) & matrix_file == 0
        return
    end
    if ~iscell(matrix_file)
        matrix_file = cellstr(matrix_file);
    end
    MatrixAll = [];
    for Nfile = 1:length(matrix_file)
        matrix_file2 = matrix_file{Nfile};
        [~,~,~,matrix_fileS] = lab_filename(matrix_file2);
        [matrix,header] = lab_read_data(fullfile(matrix_filepath,matrix_file2));
        if isempty(MatrixAll) | (size(MatrixAll,1) == size(matrix,1) & size(MatrixAll,2) == size(matrix,2))
            MatrixAll = cat(3,MatrixAll,matrix);
        end
    end
    if size(MatrixAll,3) > 1
        matrix = mean(MatrixAll,3);
        matrix_fileS = 'Average_matrix';
    end
    clearvars MatrixAll Nfile matrix_file
end

if isstruct(matrix)
    if isfield(matrix,'data')
        MatrixAll = matrix.data;
        tmp = mean(matrix.data,3);
        matrix = tmp;
        clearvars tmp
        MatrixNr = 'Mean';
    else
        return
    end
elseif size(matrix,3) > 1
    MatrixAll = matrix;
    matrix = mean(matrix,3);
    MatrixNr = 'Mean';
else
    MatrixAll = matrix;
    MatrixNr = '1';
end

if isempty(matrix) | size(matrix,1) ~= size(matrix,2)
    return
end

% set diagonal to NaN if 0 or 1
tmp = diag(matrix);
if (max(tmp) == 0 & min(tmp) == 0) | (max(tmp) == 1 & min(tmp) == 1)
    matrix(1:size(matrix,1)+1:end) = NaN;
end

if doinvisible == true
    fig1 = figure('Color',[1 1 1],'MenuBar','none','NumberTitle','off','Name','Matrix','Visible','off');
else
    fig1 = figure('Color',[1 1 1],'MenuBar','none','NumberTitle','off','Name','Matrix');
end
colormap('gray');
cmap = colormap;
colormap(flipud(cmap));
T = imagesc(matrix);
colorbar;
set(T,'UserData',matrix);

if nomenu == false
    m1 = uimenu(fig1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Save Matrix','Callback',@(hObj,evd)savematrix(T));
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig1,'Label','View');
    uimenu(m2,'Label','Clear Diagonal','Callback',@(hObj,evd)setdiagonalview(T));
    uimenu(m2,'Label','Set Threshold','Callback',@(hObj,evd)setthreshold(T));
    m21 = uimenu(m2,'Label','Histogram');
    uimenu(m21,'Label','with zeros','Callback',@(hObj,evd)histogram1(T));
    uimenu(m21,'Label','without zeros','Callback',@(hObj,evd)histogram2(T));
    m22 = uimenu(m2,'Label','Distribution');
    uimenu(m22,'Label','with zeros','Callback',@(hObj,evd)distribution1(T));
    uimenu(m22,'Label','without zeros','Callback',@(hObj,evd)distribution2(T));
    uimenu(m2,'Label','Value-Plot','Callback',@(hObj,evd)valueplot(T));
    m3 = uimenu(fig1,'Label','Edit');
    uimenu(m3,'Label','Edit Matrix','Callback',@(hObj,evd)editmatrix(T));
    uimenu(m3,'Label','Set Diagonal to NaN','Callback',@(hObj,evd)setdiagonal(T));
    uimenu(m3,'Label','Normalize','Callback',@(hObj,evd)donormalize(T));
    uimenu(m3,'Label','Make symmetrical','Callback',@(hObj,evd)dosymmetrical(T));
    uimenu(m3,'Label','Rank matrix','Callback',@(hObj,evd)doranking(T));
    uimenu(m3,'Label','Correct neighbours','Callback',@(hObj,evd)correctneighbours(T));
    uimenu(m3,'Label','Reduce to Mappings','Callback',@(hObj,evd)reduce2mappings(T));
    uimenu(m3,'Label','Shuffle Matrix','Callback',@(hObj,evd)shufflematrix(T));
    uimenu(m3,'Label','Set Min/Max','Callback',@(hObj,evd)setrange(T));
    uimenu(m3,'Label','Matrix2Binary','Callback',@(hObj,evd)matrix2binary(T));
    m4 = uimenu(fig1,'Label','Graph analysis');
    uimenu(m4,'Label','Graph analysis','Callback',@(hObj,evd)do_graphanalysis(T));
    uimenu(m4,'Label','MST','Callback',@(hObj,evd)do_mst(T));
    uimenu(m4,'Label',' show mean results','Callback',@(hObj,evd)show_graph_mean(T));
    uimenu(m4,'Label',' show single results','Callback',@(hObj,evd)show_graph_single(T));
    m5 = uimenu(fig1,'Label','Plot');
    uimenu(m5,'Label','Signal space','Callback',@(hObj,evd)plot_signalspace(T));
    uimenu(m5,'Label','Source space','Callback',@(hObj,evd)plot_sourcespace(T));
    m6 = uimenu(fig1,'Label','Select');
    uimenu(m6,'Label','Dialog','Callback',@(hObj,evd)select_matrix(T));
    uimenu(m6,'Label','Next matrix','Callback',@(hObj,evd)select_next(T));
    uimenu(m6,'Label','Previous matrix','Callback',@(hObj,evd)select_previous(T));
    uimenu(m6,'Label','Mean matrix','Callback',@(hObj,evd)select_mean(T));
end
%if exist('matrix_fileS','var')
%     lab_print_figure(fullfile(matrix_filepath,[matrix_fileS '.tif']),fig1);
%end

% helper function start
    function savematrix(T)
        matrixtmp = get(T,'UserData');
        [filename,filepath]=uiputfile('*.txt','Select File/Folder');
        if filename~=0
            lab_write_matrix(fullfile(filepath,filename),matrixtmp);
        end
    end

    function setthreshold(T)
        matrixtmp = get(T,'UserData');
        minval = min(matrixtmp(:));
        maxval = max(matrixtmp(:));
        answer = inputdlg({'Lower Threshold';'Higher Threshold'},'Threshold',[1 40;1 40],{num2str(minval);num2str(maxval)});
        if ~isempty(answer)
            minval = str2num(answer{1,1}); %#ok<ST2NM>
            maxval = str2num(answer{2,1}); %#ok<ST2NM>
            matrixtmp(matrixtmp<minval) = minval;
            matrixtmp(matrixtmp>maxval) = maxval;
            set(T,'CData',matrixtmp);
        end
    end

    function setdiagonalview(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        set(T,'CData',matrixtmp);
    end

    function editmatrix(T)
        matrixtmp = get(T,'UserData');
        matrixtmp = lab_table_dialog(matrixtmp);
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function setdiagonal(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function donormalize(T)
        matrixtmp = get(T,'UserData');
        matrixtmp = lab_normalize_matrix(matrixtmp,true);
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function dosymmetrical(T)
        matrixtmp = get(T,'UserData');
        matrixtmp = (matrixtmp + tril(matrixtmp)' + triu(matrixtmp)') ./ 2;
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function doranking(T)
        matrixtmp = get(T,'UserData');
        rankorder = inputdlg('Rank order','Rank order',[1 15],{'5'});
        rankorder = str2num(rankorder{1}); %#ok<ST2NM>
        matrixtmp = lab_rank_matrix(matrixtmp,rankorder);
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function correctneighbours(T)
        matrixtmp = get(T,'UserData');
        settings = lab_get_DIST([],true);
        if isfield(settings,'matrix') & ~isempty(settings.matrix) & ...
                size(settings.matrix,1) == size(matrixtmp,1) & ...
                size(settings.matrix,2) == size(matrixtmp,2)
            switch settings.mode
                case '1'
                    matrixtmp(settings.matrix <= settings.threshold) = 1;
                case 'max'
                    matrixtmp(settings.matrix <= settings.threshold) = max(matrixtmp(:));
            end
            set(T,'CData',matrixtmp,'Userdata',matrixtmp);
        end
    end

    function reduce2mappings(T)
        matrixtmp = get(T,'UserData');
        cfg.EXTRA.numdatachans = size(matrixtmp,1);
        Mappings = lab_load_mappings([],cfg);
        if ~isempty(Mappings) & isfield(Mappings,'mappingsChannels') & ...
                Mappings.mappingsChannels == size(matrixtmp,1);
            settings.MATRIX.Mappings = Mappings;
            settings.MATRIX.MappingsMode = questdlg('Mode','Select mode','Betweenness','Degree','Average','Average');
            R = lab_reduce_matrix2mappings(matrixtmp,settings,true);
            matrixtmp = R.matrix;
            clearvars settings R
            set(T,'CData',matrixtmp,'Userdata',matrixtmp);
            set(gca,'YLim',[0.5 size(matrixtmp,1)+0.5],'XLim',[0.5 size(matrixtmp,1)+0.5]);
        end
    end

    function shufflematrix(T)
        matrixtmp = get(T,'UserData');
        if length(unique(matrixtmp(:))) == 2 & min(matrixtmp(:)) == 0 & max(matrixtmp(:)) == 1
            numchans = size(matrixtmp,1);
            Dtmp = diag(matrixtmp);
            matrixtmp(1:numchans+1:end) = 0;
            matrixtmp = randmio_und(matrixtmp,5);
            matrixtmp(1:numchans+1:end) = Dtmp;
        else
            matrixtmp = lab_rand_matrix_fixed(matrixtmp,5);
        end
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function matrix2binary(T)
        matrixtmp = get(T,'UserData');
        if length(unique(matrixtmp(:))) == 2 & min(matrixtmp(:)) == 0 & max(matrixtmp(:)) == 1
            return
        else
            matrixtmp = lab_matrix2binary(matrixtmp);
        end
        set(T,'CData',matrixtmp,'Userdata',matrixtmp);
    end

    function do_graphanalysis(T)
        matrixtmp = get(T,'UserData');
        Graph = lab_graphanalysis(matrixtmp,header);
    end

    function do_mst(T)
        matrixtmp = get(T,'UserData');
        matrixtmp = lab_MST(matrixtmp);
        set(T,'Cdata',matrixtmp);
        [filename,filepath]=uiputfile('*.txt','Select Matrix-file to store');
        if ~isnumeric(filename);
           lab_write_matrix(fullfile(filepath,filename),matrixtmp);
        end
    end

    function show_graph_single(T)
        if isempty(Graph)
            matrixtmp = get(T,'UserData');
            Graph = lab_graphanalysis(matrixtmp,header);
        end
        Rsingle = cat(2,Graph.results,num2cell(permute(Graph.data,[3 1 2])));
        lab_table_dialog(Rsingle,[],'Graph anaylsis - Single',0);
    end

    function show_graph_mean(T)
        if isempty(Graph)
            matrixtmp = get(T,'UserData');
            Graph = lab_graphanalysis(matrixtmp,header);
        end
        Rmean = cat(2,Graph.resultsmean,num2cell(Graph.datamean));
        lab_table_dialog(Rmean,[],'Graph anaylsis - Mean',0);
    end

    function histogram1(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        settings.excludezeros = false;
        settings.combine = true;
        lab_plot_histogram(matrixtmp(:),settings);
    end

    function histogram2(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        settings.excludezeros = true;
        settings.combine = true;
        lab_plot_histogram(matrixtmp(:),settings);
    end

    function distribution1(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        lab_plot_distribution(matrixtmp(:),true);
    end

    function distribution2(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        matrixtmp(matrixtmp==0) = NaN;
        lab_plot_distribution(matrixtmp(:),true);
    end

    function valueplot(T)
        matrixtmp = get(T,'UserData');
        matrixtmp(1:size(matrixtmp,1)+1:end) = NaN;
        matrixtmp(matrixtmp==0) = NaN;
        f = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none','Name','Sorted Values - Plot');
        Tm1 = uimenu(f,'Label','File');
        uimenu(Tm1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(Tm1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(Tm1,'Label','Close','Callback','close;');
        plot(sort(matrixtmp(:)));
    end

    function setrange(T)
        A = get(T,'Parent');
        CLim = get(A,'CLim');
        settings.minval = CLim(1);
        settings.maxval = CLim(2);
        Prompt = {'Min value','minval';'Max value','maxval'};
        Formats.type = 'edit';
        Formats.format = 'float';
        Formats.size = 30;
        Formats = [Formats Formats];
        [settings,Cancelled] = inputsdlg(Prompt,'Range',Formats,settings);
        if Cancelled ~= 1
            CLim = [settings.minval settings.maxval];
            set(A,'CLim',CLim);
        end
    end

    function plot_signalspace(T)
        matrixtmp = get(T,'UserData');
        lab_plot_elec(matrixtmp);
    end

    function plot_sourcespace(T)
        matrixtmp = get(T,'UserData');
        lab_plot_IS(matrixtmp);
    end
    
    function select_matrix(T)
        if size(MatrixAll,3) > 1
            strlist = cellstr(num2str((1:size(MatrixAll,3))'));
            strlist = cat(1,cellstr('Mean'),strlist);
            defaultstr = find(strcmp(strlist,MatrixNr));
            selection = listdlg('PromptString','Matrix:','SelectionMode','single', ...
                'ListString',strlist,'InitialValue',defaultstr,'ListSize',[80 200]);
            if ~isempty(selection)
                if selection == 1
                    set(T,'CData',mean(MatrixAll,3));
                    MatrixNr = 'Mean';
                else
                    set(T,'CData',MatrixAll(:,:,selection-1));
                    MatrixNr = num2str(selection-1);
                end
            end
        end
    end
    
    function select_next(T)
        if size(MatrixAll,3) > 1
            Nr = str2num(MatrixNr); %#ok<ST2NM>
            if isempty(Nr)
                Nr = 1;
            end
            if Nr < size(MatrixAll,3)
                Nr = Nr + 1;
            end
            set(T,'CData',MatrixAll(:,:,Nr));
            MatrixNr = num2str(Nr);
        end
    end
    
    function select_previous(T)
        if size(MatrixAll,3) > 1
            Nr = str2num(MatrixNr); %#ok<ST2NM>
            if isempty(Nr)
                Nr = 1;
            end
            if Nr > 1
                Nr = Nr - 1;
            end
            set(T,'CData',MatrixAll(:,:,Nr));
            MatrixNr = num2str(Nr);
        end
    end
    
    function select_mean(T)
        if size(MatrixAll,3) > 1
            set(T,'CData',mean(MatrixAll,3));
            MatrixNr = 'Mean';
        end
    end
% helper function end

end