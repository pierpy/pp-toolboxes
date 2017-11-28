% For xls-file with results and outcomes in last row, this function will
% split the file in one xls-file per outcome
%
% lab_split_by_outcome
%
% written by F. Hatz

function lab_split_by_outcome

% Read xls-file
[headertmp,~,~,~,cfg] = lab_read_statistics([],1,0,0,0,0);

if ~isempty(headertmp)
    headertmp = headertmp';
    header.file = cfg.file;
    header.path = cfg.path;
    try
        tmp = headertmp(2:end,2:end);
        tmp2 = find(cellfun(@isempty,tmp));
        if ~isempty(tmp2)
            for i = 1:length(tmp2);
                tmp{tmp2(i)} = NaN;
            end
        end
        data = cell2mat(tmp);
    catch %#ok<CTCH>
        disp('Input file has invalid / non numeric data')
        return
    end
else
    return
end
if isempty(data)
    return
end

if ~isfield(cfg,'numresults')
    cfg.numresults = 1;
end

strlist = headertmp(1,end-cfg.numresults+1:end);
selection = listdlg('PromptString','Select outcome:','SelectionMode','single','ListString',strlist);
result = data(:,selection + end-cfg.numresults);

iresult = unique(result);

for i = 1:length(iresult)
    tmp = find(result == iresult(i));
    tmp = [1;tmp+1];
    xlsout = headertmp(tmp,:);
    if size(xlsout,2) > 255
        fileout = fullfile(header.path,[header.file '_' num2str(i) '.xlsx']);
    else
        fileout = fullfile(header.path,[header.file '_' num2str(i) '.xls']);
    end
    lab_write_xls(fileout,xlsout');
    clearvars xlsout tmp
end

return