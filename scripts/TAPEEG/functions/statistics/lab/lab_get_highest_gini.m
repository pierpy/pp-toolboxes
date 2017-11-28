function lab_get_highest_gini

[tmp1,tmp2] = uigetfile('*.xls;*.xlsx','Select File with Gini importances');
if isnumeric(tmp1)
    return
end
Gini_file = fullfile(tmp2,tmp1);
clearvars tmp1 tmp2
Gini = lab_read_xls(Gini_file);

[tmp1,tmp2] = uigetfile('*.xls;*.xlsx','Select Input data');
if isnumeric(tmp1)
    return
end
Data_file = fullfile(tmp2,tmp1);
clearvars tmp1 tmp2
inputdata = lab_read_xls(Data_file);

settings.number = 10;
Prompt = {'Number of variables to include','number'};
Formats.type = 'edit';
[settings,Cancelled] = inputsdlg(Prompt,'Set number',Formats,settings);
if Cancelled == 1
    return
end

[~,Isort] = sort(Gini,'descend');
if length(Isort) > settings.number
    Isort = Isort(1:settings.number);
end

dataout = inputdata(1,:);
dataout = cat(1,dataout,inputdata(Isort+1,:));
dataout = cat(1,dataout,inputdata(end,:));

[~,Data_filepath,~,Data_fileS] = lab_filename(Data_file);
lab_write_xls(fullfile(Data_filepath,[Data_fileS '_selected.xlsx']),dataout);

end

