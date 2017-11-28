function lab_write_struct2xls(Filename,Var)

Names = fieldnames(Var);
Names = sort(Names);
List = {};
List2 = {};
for i = 1:length(Names)
    if size(Var.(Names{i}),1) == 1 | size(Var.(Names{i}),2) == 1
        tmp = Var.(Names{i});
        if isnumeric(tmp)
            tmp = num2cell(tmp(:)');
        elseif ischar(tmp)
            tmp = cellstr(tmp);
        elseif ~iscell(tmp)
            tmp = cellstr('invalid data');
        end
        List{end+1,1} = Names{i}; %#ok<AGROW>
        if ~isempty(tmp)
            List(end,2:length(tmp)+1) = tmp;
        end
    elseif isempty(Var.(Names{i}))
        List{end+1,1} = Names{i}; %#ok<AGROW>
    else
        tmp = Var.(Names{i});
        if isnumeric(tmp)
            tmp = num2cell(tmp);
        elseif ischar(tmp)
            tmp = cellstr(tmp);
        elseif ~iscell(tmp)
            tmp = cellstr('invalid data');
        end
        if ~isempty(tmp)
            List2{1,end+1} = tmp; %#ok<AGROW>
            List2{2,end} = Names{i};
        end
    end
end

[~,Filepath,~,FilenameS] = lab_filename(Filename);
warning off %#ok<WNOFF>
if exist(fullfile(Filepath,[FilenameS '_bad.xlsx']),'file')
    delete(fullfile(Filepath,[FilenameS '_bad.xlsx']));
end
lab_write_xls(fullfile(Filepath,[FilenameS '_bad.xlsx']),List,'BAD');
for i = 1:size(List2,2)
    lab_write_xls(fullfile(Filepath,[FilenameS '_bad.xlsx']),List2{1,i},List2{2,i});
end
warning on %#ok<WNON>