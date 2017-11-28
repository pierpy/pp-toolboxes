% Scatter-Plot Connectivity values against distance
% Input is matrix-file and electrodes-locations
%
% written by F. Hatz 2014

function lab_plot_matrix2distance

settings.Matrix = [];
settings.Matrix_file = '';
settings.locs = [];
settings.Marker.size = 36;
settings.Marker.marker = 'o';
settings.statistics = 'None';

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Matrix-File','Matrix_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.txt','Matrix-file'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = 300;
Formats(end,1).callback = {@get_matrix,{'Matrix','Matrix_file'},'Matrix_file'};

Prompt(end+1,:) = {'Matrix','Matrix'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_matrix,{'Matrix','Matrix_file'},'Matrix','Matrix_file'};

Prompt(end+1,:) = {'LOCS-file','locs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_locs,'locs','locs','Matrix'};

Prompt(end+1,:) = {'Marker','Marker'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@get_Marker,'Marker','Marker'};

Prompt(end+1,:) = {'Statistics','statistics'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'None','Pearson','Kendall','Spearman'};

[settings,Cancelled] = inputsdlg(Prompt,'Plot Matrix & Distance',Formats,settings);
if isempty(settings) | Cancelled == 1
    return
end
pause(0.2);

if isempty(settings.Matrix) | isempty(settings.locs)
    return
end

matrix = lab_rm_diagonal(settings.Matrix);
if ~isempty(settings.Matrix_file)
    [~,~,~,Title] = lab_filename(settings.Matrix_file);
    Title = regexprep(Title,'_',' ');
else
    Title = 'Matrix & Distance';
end

distance = lab_distance(settings.locs);
distance = lab_rm_diagonal(distance);

if strcmp(settings.statistics,'Kendall') | strcmp(settings.statistics,'Spearman') | strcmp(settings.statistics,'Pearson')
    disp(['    calculate statistics: ' settings.statistics])
    [coeff,pval]=corr(distance(:),matrix(:),'type',settings.statistics);
else
    pval = [];
    coeff = [];
end

if size(matrix,2) ~= size(distance,2)
    disp('Abort: Matrix and Electrode-file not matching')
    return
end

figure1 = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none','Name',['Matrix & Distance : ' Title]);
m1 = uimenu(figure1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
        
scatter(distance(:),matrix(:),settings.Marker.size,[0 0 0],settings.Marker.marker);
if ~isempty(pval)
    title ([settings.statistics ' coefficient: ' num2str(coeff,'%0.3f') ' p-value: ' num2str(pval,'%0.3f')]);
end

end

function [Matrix,Matrix_file] = get_matrix(Matrix_file)
    if exist(Matrix_file,'file')
        tmp = lab_read_data(Matrix_file);
        if ~isempty(tmp) & isnumeric(tmp) & size(tmp,1) == size(tmp,2)
            Matrix = tmp;
        end
    else
        Matrix = [];
        Matrix_file = '';
    end
end

function [Matrix,Matrix_file] = load_matrix(Matrix,Matrix_file)
    Matrix_fileB = Matrix_file;
    MatrixB = Matrix;
    [Matrix,Matrix_file] = lab_load_matrix(Matrix);
    if corr(Matrix(:),MatrixB(:)) == 1 & isempty(Matrix_file)
        Matrix_file = Matrix_fileB;
    end
end

function locs = load_locs(locs,Matrix)
    if ~isempty(Matrix)
        numchans = size(Matrix,1);
        locs = lab_load_locs(locs,[],numchans);
    else
        locs = lab_load_locs(locs);
    end
end

function Marker = get_Marker(Marker)

if isempty(Marker)
    Marker.size = 36;
    Marker.marker = 'o';
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Marker','marker'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'o';'.';'*';'x'};

Prompt(end+1,:) = {'',[]};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Size','size'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).size = 50;
Formats(end,1).limits = [0 999];

[Marker,Cancelled] = inputsdlg(Prompt,'Marker',Formats,Marker);
if Cancelled == 1
    Marker = [];
end

end