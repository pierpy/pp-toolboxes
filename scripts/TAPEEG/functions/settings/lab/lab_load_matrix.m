% Load Matrix info
%
%
% cfg    = structure with config (optional)
%
% written by F. Hatz 2012

function [Matrix,Matrix_file] = lab_load_matrix(Matrix,numchans,doMST)

Diagonal = [];
if ~exist('doMST','var')
    doMST = false;
end
if ~exist('numchans','var')
    numchans = [];
end
if exist('Matrix','var') & ~isempty(Matrix) & ischar(Matrix) & exist(Matrix,'file')
    tmp = lab_read_data(Matrix);
    if ~isempty(tmp) & isnumeric(tmp) & size(tmp,1) == size(tmp,2)
        settings.Matrix = tmp;
        settings.Matrix_file = Matrix;
        numchans = size(settings.Matrix,1);
    else
        settings.Matrix = [];
        settings.Matrix_file = '';
    end
elseif exist('Matrix','var') & ~isempty(Matrix) & isnumeric(Matrix) & size(Matrix,1) == size(Matrix,2)
    Diagonal = diag(Matrix);
    if (max(Diagonal) == 0 & min(Diagonal) == 0) | (max(Diagonal) == 1 & min(Diagonal) == 1)
        Matrix(1:size(Matrix,1)+1:end) = NaN;
    end
    settings.Matrix = Matrix;
    numchans = size(Matrix,1);
    settings.Matrix_file = '';
elseif ~isempty(numchans) & isnumeric(numchans) & numchans > 1
    settings.Matrix = zeros(numchans,numchans);
    settings.Matrix_file = '';
else
    settings.Matrix = [];
    settings.Matrix_file = '';
end
settings.nchans = numchans;

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Matrix-File','Matrix_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.txt','TXT-file(s)';'*.xls;*.xlsx','Excel-File(s)'};
Formats(end,1).limits = [0 5];
Formats(end,1).callback =  {@read_matrix,{'Matrix','Matrix_file','nchans'},'Matrix_file'};
Formats(end,1).size = [360 40];
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Number channels','nchans'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).size = 40;
Formats(end,1).callback = {@set_numberchannels,'Matrix','Matrix','nchans'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Create',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@create_matrix,{'Matrix','nchans'},'nchans'};

Prompt(end+1,:) = {'','Matrix'};
Formats(end+1,1).type = 'image';
Formats(end,1).format = 'matrix';
Formats(end,1).size = [300 300];
Formats(end,1).items = {repmat((1:-1/255:0)',1,3)};
Formats(end,1).enable = 'inactive';
Formats(end,1).callback = {@set_point,'Matrix','Matrix','@Matrix'};
Formats(end,1).span = [5 1];

Prompt(end+1,:) = {'Edit',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@lab_table_dialog,'Matrix','Matrix'};

Prompt(end+1,:) = {'Shuffle',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@shuffle_matrix,'Matrix','Matrix'};

if doMST == false
    Prompt(end+1,:) = {'MST',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@convert_MST,'Matrix','Matrix'};
end

Prompt(end+1,:) = {'Save',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@save_matrix,'Matrix_file','Matrix','Matrix_file'};

[settings,Cancelled] = inputsdlg(Prompt,'Load Matrix',Formats,settings);
if Cancelled == 1
    Matrix = [];
    Matrix_file = '';
    return
end
Matrix = settings.Matrix;
if ischar(settings.Matrix_file)
    Matrix_file = settings.Matrix_file;
elseif iscell(settings.Matrix_file) & ~isempty(settings.Matrix_file)
    Matrix_file = settings.Matrix_file{1};
else
    Matrix_file = '';
end
if ~isempty(Diagonal) & length(Diagonal) == size(Matrix,1)
    Matrix(1:size(Matrix,1)+1:end) = Diagonal;
end

    function [Matrix,Matrix_file,Nchans] = read_matrix(Matrix_file)
        if isempty(Matrix_file)
            Matrix = [];
            return
        end
        if ischar(Matrix_file)
            Matrix_file = cellstr(Matrix_file);
        end
        MatrixAll = [];
        for Nfile = 1:length(Matrix_file)
            Matrix_file2 = Matrix_file{Nfile};
            Matrix = lab_read_data(Matrix_file2);
            if isempty(MatrixAll) | (size(MatrixAll,1) == size(Matrix,1) & size(MatrixAll,2) == size(Matrix,2))
                MatrixAll = cat(3,MatrixAll,Matrix);
            end
        end
        if size(MatrixAll,3) > 1
            Matrix = mean(MatrixAll,3);
            [~,Pathtmp] = lab_filename(Matrix_file{end});
            Matrix_file = {fullfile(Pathtmp,'Average_matrix.txt')};
        else
            Matrix_file = Matrix_file{1};
        end
        clearvars MatrixAll Nfile
        if ~isnumeric(Matrix) | size(Matrix,1) ~= size(Matrix,2)
            Matrix = [];
            Matrix_file = cellstr('');
            Nchans = [];
            return
        end
        if doMST == true
            Matrix = lab_MST(Matrix);
        end
        Diagonal = diag(Matrix);
        if (max(Diagonal) == 0 & min(Diagonal) == 0) | (max(Diagonal) == 1 & min(Diagonal) == 1)
            Matrix(1:size(Matrix,1)+1:end) = NaN;
        end
        Nchans = size(Matrix,1);
    end

    function [Matrix,Nchans] = create_matrix(Nchans)
        Matrix = lab_generate_matrix(Nchans);
    end

    function Matrix = shuffle_matrix(Matrix)
        if length(unique(Matrix(:))) == 2 & min(Matrix(:)) == 0 & max(Matrix(:)) == 1
            numchans = size(Matrix,1);
            Dtmp = diag(Matrix);
            Matrix(1:numchans+1:end) = 0;
            Matrix = randmio_und(Matrix,5);
            Matrix(1:numchans+1:end) = Dtmp;
        else
            Matrix = lab_rand_matrix_fixed(Matrix,5);
        end
    end

    function Matrix = convert_MST(Matrix)
        if ~isempty(Diagonal) & length(Diagonal) == size(Matrix,1)
            Matrix(1:size(Matrix,1)+1:end) = Diagonal;
        end
        Matrix = lab_MST(Matrix);
    end

    function Matrix_file = save_matrix(Matrix,Matrix_file)
        [filename,filepath] = uiputfile('*.txt','Save to');
        if filename ~= 0
            if ~isempty(Diagonal)
                Matrix(1:size(Matrix,1)+1:end) = Diagonal;
            end
            dlmwrite(fullfile(filepath,filename),Matrix,'delimiter','\t','precision', 6);
            Matrix_file = fullfile(filepath,filename);
        end
    end

    function Matrix = set_point(Matrix,Mhandle)
        Mtmp = tril(Matrix) - triu(Matrix)';
        if max(abs(Mtmp(:))) < 10^-6
            symmetrical = true;
        else
            symmetrical = false;
        end
        T = get(Mhandle,'CurrentPoint');
        Tx = ceil(T(1,1)/300 * size(Matrix,1));
        Ty = ceil(T(1,2)/300 * size(Matrix,2));
        Munique = unique(Matrix(:));
        Munique = Munique(~isnan(Munique));
        if length(Munique) == 2 & min(Munique) == 0 & max(Munique) == 1
            if Matrix(Ty,Tx) == 0
                Matrix(Ty,Tx) = 1;
                if symmetrical == true
                    Matrix(Tx,Ty) = 1;
                end
            else
                Matrix(Ty,Tx) = 0;
                if symmetrical == true
                    Matrix(Tx,Ty) = 0;
                end
            end
        elseif length(Munique) == 1 & max(Munique) == 0
            Matrix(Ty,Tx) = 1;
            if symmetrical == true
                Matrix(Tx,Ty) = 1;
            end
        else
            Mtmp = inputdlg('Set value','Value',[1 40],cellstr(num2str(Matrix(Ty,Tx))));
            Matrix(Ty,Tx) = str2num(Mtmp{1});
            if symmetrical == true
                Matrix(Tx,Ty) = Matrix(Ty,Tx);
            end
        end
    end

    function Matrix = set_numberchannels(Matrix,nchans)
        if nchans <= size(Matrix,1)
            Matrix = Matrix(1:nchans,1:nchans);
        else
            Matrix(nchans,nchans) = min(Matrix(:));
        end
    end

end