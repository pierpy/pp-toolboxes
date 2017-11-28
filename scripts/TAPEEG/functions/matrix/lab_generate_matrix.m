function [matrix,numchans] = lab_generate_matrix(numchans)

if ~exist('numchans','var')
    settings.numchans = [];
else
    settings.numchans = numchans;
end

settings.method = 'shift';
Formats = {};
Prompt = cell(0,2);

Prompt(end+1,:) = {'Number of channels', 'numchans'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 99999];
Formats(end,1).size = 60;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Method','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
if isempty(numchans) | numchans == 78
    Formats(end,1).items = {'shift','modules','gong binary','gong weighted','zeros'};
else
    Formats(end,1).items = {'shift','modules'};
end

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Format','format'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'full','lower tirangle','upper triangle'};

[settings,Cancelled] = inputsdlg(Prompt,'Matrix Template',Formats,settings);
if Cancelled == 1
    matrix = zeros(numchans,numchans);
    return
elseif isfield(settings,'numchans')
    numchans = settings.numchans;
end
matrix = zeros(numchans,numchans);
switch settings.method
    case 'shift'
        tmp = inputdlg({'weigths';'step'},'Shift matrix',[1 20;1 20],{num2str(1);num2str(1)});
        weight = str2num(tmp{1,1});
        step = str2num(tmp{2,1});
        matrix(2:(step*(numchans+1)):end) =weight;
    case 'modules'
        tmp = inputdlg({'number of modules','weights'},'Modules matrix',[1 20;1 20],{num2str(3),num2str(1)});
        modules = round(str2num(tmp{1,1}));
        weight = str2num(tmp{2,1});
        nmodules = ceil(numchans / modules);
        for i = 1:modules-1
            if i == 1
                tmp = 1;
            else
                tmp = (i-1)*nmodules;
            end
            matrix((i-1)*nmodules+1:i*nmodules,tmp) =weight;
        end
        matrix(i*nmodules+1:end,i*nmodules) =weight;
    case 'gong binary'
        matrix = lab_gongmatrix;
    case 'gong weighted'
        matrix = lab_gongmatrixW;
    case 'zeros'
        matrix = zeros(numchans,numchans);
end
if settings.format == 2
    matrix = tril(matrix);
elseif settings.format == 3
    matrix = triu(matrix);
else
    tmp = matrix(1:size(matrix,1)+1:end);
    matrix = tril(matrix);
    matrix = matrix + matrix';
    matrix(1:size(matrix,1)+1:end) = tmp;
end