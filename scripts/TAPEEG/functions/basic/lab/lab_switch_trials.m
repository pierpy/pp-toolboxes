function lab_switch_trials
    
[Filename,Filepath] = uigetfile('*.xls;*.xlsx','Select xls-file');
if isnumeric(Filename)
    return
end

Data = lab_read_xls(fullfile(Filepath,Filename));
if isempty(Data)
    return
end

title = Data(1,2:end);
for i = 1:length(title)
    tmp = title{i};
    tmp2 = strfind(tmp,'_');
    if isempty(tmp2)
        tmp2 = strfind(tmp,' ');
    end
    if ~isempty(tmp2)
        tmp = [tmp(tmp2(end)+1:end) '_' tmp(1:tmp2(end)-1)];
        title{i} = tmp;
    end
end
Data(1,2:end) = title;
[~,~,~,FilenameS] = lab_filename(Filename);
lab_write_xls(fullfile(Filepath,[FilenameS '_trsp.xlsx']),Data);

end
