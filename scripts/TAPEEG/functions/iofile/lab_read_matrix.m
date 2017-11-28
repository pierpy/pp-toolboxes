function [matrix,header,cfg] = lab_read_matrix(Filename,cfg,doaverage,forceselection)

if ~exist('forceselection','var')
    forceselection = false;
end
if ~exist('doaverage','var')
    doaverage = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('Filename','var') | isempty(Filename)
    [Matrix_file,Matrix_filepath]=uigetfile('*.*','Select Matrix file (Text/Excel/BrainWave/Mat-Container)','MultiSelect','on');
    if isnumeric(Matrix_file)
        matrix = [];
        header = [];
        return
    end
    if iscell(Matrix_file)
        for i = 1:length(Matrix_file)
            Filename{i} = fullfile(Matrix_filepath,Matrix_file{i});
        end
    else
        Filename = fullfile(Matrix_filepath,Matrix_file);
    end
    clearvars Matrix_file Matrix_filepath
end
if ~iscell(Filename)
    Filename = cellstr(Filename);
end
[~,~,Matrix_format] = lab_filename(Filename{1});

MatrixAll = [];
header = [];
for Nfile = 1:length(Filename)
    [Matrix_file,Matrix_filepath] = lab_filename(Filename{Nfile});
    if strcmp(Matrix_format,'xls') | strcmp(Matrix_format,'xlsx')
        matrix = lab_read_xls(fullfile(Matrix_filepath,Matrix_file));
        if isnumeric(matrix)
            matrix = num2cell(matrix);
        end
        if isempty(matrix)
            return
        end
        if ischar(matrix{1,1})
            matrix = matrix(2:end,2:end);
        end
        try
            matrix = cell2mat(matrix);
        catch %#ok<CTCH>
            return
        end
        [patient,cfg] = lab_subjectname(fullfile(Matrix_filepath,Matrix_file),cfg);
        if ~isfield(header,'subjects')
            header.subjects = cellstr(patient);
        else
            header.subjects = cat(2,header.subjects,cellstr(patient));
        end
    elseif strcmp(Matrix_format,'mat')
        MAT = load(fullfile(Matrix_filepath,Matrix_file));
        if isfield(MAT,'result') & ~isempty(MAT.result)
            if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'matrixvar')
                cfg.EXTRA.matrixvar = '';
                forceselection = true;
            end
            if forceselection == true
                Fnames = fieldnames(MAT.result);
                Idx = [];
                for i = 1:length(Fnames)
                    if size(MAT.result.(Fnames{i}),1) == size(MAT.result.(Fnames{i}),2) & size(MAT.result.(Fnames{i}),1) > 1
                        Idx = [Idx i]; %#ok<AGROW>
                    end
                end
                if ~isempty(Idx)
                    Fnames = Fnames(Idx);
                    strdefault = find(strcmp(Fnames,cfg.EXTRA.matrixvar));
                    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
                        selection = listdlg('PromptString','Select Matrix:','SelectionMode','single', ...
                            'ListString',Fnames,'InitialValue',strdefault,'ListSize',[100 100]);
                        if ~isempty(selection)
                            Mvar = Fnames{selection};
                            cfg.EXTRA.matrixvar = Mvar;
                        else
                            Mvar = '';
                        end
                    else
                        Mvar = Fnames{1};
                        cfg.EXTRA.matrixvar = Mvar;
                    end
                else
                    Mvar = '';
                end
                clearvars Fnames selection Idx i
            else
                Mvar = cfg.EXTRA.matrixvar;
            end
            if ~isempty(Mvar) & isfield(MAT.result,Mvar)
                matrix = MAT.result.(Mvar);
                if isfield(MAT.result,[Mvar '_timestamp'])
                    if ~isfield(header,'timestamp')
                        header.timestamp = MAT.result.([Mvar '_timestamp'])(:)';
                    else
                        header.timestamp = cat(2,header.timestamp,MAT.result.([Mvar '_timestamp'])(:)');
                    end
                end
                if isfield(MAT.result,'locs')
                    header.locs = MAT.result.locs;
                end
                if isfield(MAT,'patient')
                    patient = MAT.patient;
                else
                    [patient,cfg] = lab_subjectname(fullfile(Matrix_filepath,Matrix_file),cfg);
                end
                if ~isfield(header,'subjects')
                    header.subjects = repmat(cellstr(patient),1,size(matrix,3));
                else
                    header.subjects = cat(2,header.subjects,repmat(cellstr(patient),1,size(matrix,3)));
                end
                if isfield(MAT.result,'freqband')
                    header.freqband = MAT.result.freqband;
                else
                    tmp = strfind(Filename,'Conn_F');
                    if ~isempty(tmp)
                        tmp = Filename(tmp+6:end-4);
                        tmp2 = strfind(tmp,'F');
                        if isempty(tmp2)
                            tmp2 = strfind(tmp,'_');
                        end
                        if ~isempty(tmp2)
                            tmp3 = tmp(tmp2(1)+1:end);
                            tmp4 = strfind(tmp3,'_');
                            if isempty(tmp4)
                                header.freqband = [str2num(tmp(1:tmp2(1)-1)) str2num(tmp(tmp2(1)+1:end))]; %#ok<ST2NM>
                            else
                                header.freqband = [str2num(tmp(1:tmp2(1)-1)) str2num(tmp3(1:tmp4(1)-1))]; %#ok<ST2NM>
                            end
                        end
                    end
                end
                if isfield(MAT.result,'channels')
                    header.channels = MAT.result.channels;
                end
            end
        elseif isfield(MAT,'Result') & ~isempty(MAT.Result)
            if ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'matrixvar')
                cfg.EXTRA.matrixvar = '';
                forceselection = true;
            end
            if forceselection == true
                Fnames = fieldnames(MAT.Result);
                Idx = [];
                for i = 1:length(Fnames)
                    if isstruct(MAT.Result.(Fnames{i})) & isfield(MAT.Result.(Fnames{i}),'matrix') & ~isempty(MAT.Result.(Fnames{i}).matrix)
                        Idx = [Idx i]; %#ok<AGROW>
                    end
                end
                if ~isempty(Idx)
                    Fnames = Fnames(Idx);
                    Fnames = cat(1,cellstr('ALL'),Fnames(:));
                    strdefault = find(strcmp(Fnames,cfg.EXTRA.matrixvar));
                    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
                        selection = listdlg('PromptString','Select Matrix:','SelectionMode','single', ...
                            'ListString',Fnames,'InitialValue',strdefault,'ListSize',[100 100]);
                        if ~isempty(selection)
                            Mvar = Fnames{selection};
                            cfg.EXTRA.matrixvar = Mvar;
                        else
                            Mvar = '';
                        end
                    else
                        Mvar = Fnames{1};
                        cfg.EXTRA.matrixvar = Mvar;
                    end
                else
                    Mvar = '';
                end
                clearvars Fnames selection Idx i
            else
                Mvar = cfg.EXTRA.matrixvar;
            end
            if ~isempty(Mvar) & isfield(MAT.Result,Mvar)
                if isstruct(MAT.Result.(Mvar).matrix)
                    Mnames = fieldnames(MAT.Result.(Mvar).matrix);
                    selection = listdlg('PromptString','Select Matrix:','SelectionMode','single', ...
                        'ListString',Mnames,'ListSize',[100 100]);
                    if ~isempty(selection)
                        matrix = MAT.Result.(Mvar).matrix.(Mnames{selection(1)});
                    else
                        matrix = [];
                        header = [];
                        return
                    end
                elseif isnumeric(MAT.Result.(Mvar).matrix)
                    matrix = MAT.Result.(Mvar).matrix;
                else
                    matrix = [];
                    header = [];
                    return
                end
                if isfield(MAT.Result.(Mvar),'subjects')
                    if exist('Mnames','var')
                        subjects = MAT.Result.(Mvar).subjects.(Mnames{selection(1)})(:)';
                    else
                        subjects = MAT.Result.(Mvar).subjects(:)';
                    end
                    if ~isfield(header,'subjects')
                        header.subjects = subjects;
                    else
                        header.subjects = cat(2,header.subjects,subjects);
                    end
                end
                if isfield(MAT.Result.(Mvar),'channels')
                    header.channels = MAT.Result.(Mvar).channels;
                end
            elseif ~isempty(Mvar) & strcmp(Mvar,'ALL')
                Fnames = fieldnames(MAT.Result);
                Idx = [];
                for i = 1:length(Fnames)
                    if isstruct(MAT.Result.(Fnames{i})) & isfield(MAT.Result.(Fnames{i}),'matrix') & ~isempty(MAT.Result.(Fnames{i}).matrix)
                        Idx = [Idx i]; %#ok<AGROW>
                    end
                end
                if ~isempty(Idx)
                    Fnames = Fnames(Idx);
                else
                    Fnames = [];
                end
                if ~isempty(Fnames)
                    matrix = [];
                    header.subjects = {};
                    for i = 1:length(Fnames)
                        matrix = cat(3,matrix,MAT.Result.(Fnames{i}).matrix);
                        if isfield(MAT.Result.(Fnames{i}),'subjects')
                            header.subjects = cat(2,header.subjects,MAT.Result.(Fnames{i}).subjects(:)');
                        end
                        if isfield(MAT.Result.(Fnames{i}),'channels')
                            header.channels = MAT.Result.(Fnames{i}).channels;
                        end
                    end
                end
            end
            tmp = strfind(Matrix_file,'_F');
            if ~isempty(tmp) & ~isnan(str2double(Matrix_file(tmp+2)))
                tmp = Matrix_file(tmp+2:end-4);
                tmp2 = strfind(tmp,'F');
                if isempty(tmp2)
                    tmp2 = strfind(tmp,'_');
                end
                if ~isempty(tmp2)
                    tmp3 = tmp(tmp2(1)+1:end);
                    tmp4 = strfind(tmp3,'_');
                    if isempty(tmp4)
                        header.freqband = [str2num(tmp(1:tmp2(1)-1)) str2num(tmp(tmp2(1)+1:end))]; %#ok<ST2NM>
                    else
                        header.freqband = [str2num(tmp(1:tmp2(1)-1)) str2num(tmp3(1:tmp4(1)-1))]; %#ok<ST2NM>
                    end
                end
            elseif strcmp(Mvar(1),'F')
                tmp = strfind(Mvar,'F');
                if length(tmp) == 2
                    header.freqband = [str2num(Mvar(2:tmp(2)-1)) str2num(Mvar(tmp(2)+1:end))]; %#ok<ST2NM>
                end
            end
        else
            matrix = [];
        end
    else
        [matrix,headertmp,cfg] = lab_read_txt(fullfile(Matrix_filepath,Matrix_file),cfg);
        if isstruct(matrix) & isfield(headertmp,'subjects')
            if ~isfield(header,'subjects')
                header.subjects = headertmp.subjects;
            else
                header.subjects = cat(2,header.subjects,headertmp.subjects);
            end
        else
            [patient,cfg] = lab_subjectname(fullfile(Matrix_filepath,Matrix_file),cfg);
            if ~isfield(header,'subjects')
                header.subjects = cellstr(patient);
            else
                header.subjects = cat(2,header.subjects,cellstr(patient));
            end
        end
    end
    if ~exist('matrix','var')
        matrix = [];
        header = [];
        return
    end
    if isempty(MatrixAll) | (size(MatrixAll,1) == size(matrix,1) & size(MatrixAll,2) == size(matrix,2))
        MatrixAll = cat(3,MatrixAll,matrix);
    end
end
if isempty(MatrixAll) | ~isnumeric(MatrixAll) | ~size(MatrixAll,2) == size(MatrixAll,1)
    disp('    No valid matrix information')
    matrix = [];
    header = [];
    return
end
if size(MatrixAll,3) == 1
    cfg.EEG_file = Matrix_file;
    cfg.patient = header.subjects{1};
    matrix = MatrixAll;
elseif doaverage == false
    cfg.EEG_file = 'MultiMatrix.txt';
    cfg.patient = '';
    matrix = MatrixAll;
else
    matrix = mean(MatrixAll,3);
    cfg.EEG_file = 'Average_Matrix.txt';
    if length(unique(header.subjects)) == 1
        header.subjects = header.subjects(1);
    else
        header.subjects = cellstr('AverageSubjects');
    end
    cfg.patient = header.subjects{1};
end
cfg.EEG_filepath = Matrix_filepath;
header.datatype = 'matrix';

end