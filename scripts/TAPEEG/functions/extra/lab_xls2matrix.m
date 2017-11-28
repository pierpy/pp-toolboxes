function lab_xls2matrix

[XLS_file,XLS_filepath] = uigetfile('*.xls;*.xlsx','Select XLS-File');
if isnumeric(XLS_file)
    return
end

datainput = lab_read_xls(fullfile(XLS_filepath,XLS_file));
Mtranspose = questdlg('Are Matrix values in Columns or Rows?','Values in column/row','Cancel','Row','Column','Column');
if isempty(Mtranspose) | strcmp(Mtranspose,'Cancel')
    return
elseif strcmp(Mtranspose,'Row')
    datainput = datainput';
end
datainput = correctheader(datainput);
measures = datainput(1,2:end)';
for i = 1:length(measures)
    measures{i} = regexprep(measures{i},{'::',':',';','/','\',','},'_');
end
matrixvalues = datainput(2:end,1); %#ok<NASGU>
data = cell2mat(datainput(2:end,2:end));

settings = [];
[settings,skipprocessing] = lab_set_xls2matrix(settings,measures);
if skipprocessing == 1
    return
end
data = data(:,settings.selection);
measures = measures(settings.selection);
if strcmp(settings.invertvalues,'1-value')
    data = 1-data;
elseif strcmp(settings.invertvalues,'1/value')
    data = data.^-1;
end

numvalues = size(data,1);
if floor((numvalues*2)^0.5) == ((numvalues*2) - floor((numvalues*2)^0.5))^0.5
    Nmatrix = floor((numvalues*2)^0.5);
    hasdiagonal = true;
elseif ceil((numvalues*2)^0.5) == ((numvalues*2) - ceil((numvalues*2)^0.5))^0.5
    Nmatrix = ceil((numvalues*2)^0.5);
    hasdiagonal = false;
else
    disp('Abort: Number of matrix values are not corresponding to a matrix')
    return
end

Output_filepath = fullfile(XLS_filepath,settings.folder);
[~,~,~,Output_fileS] = lab_filename(XLS_file);
warning off %#ok<WNOFF>
mkdir(Output_filepath);
warning on %#ok<WNON>

for i = 1:size(data,2)
    matrixout = zeros(Nmatrix,Nmatrix);
    T = 1;
    if hasdiagonal == true
        for m = 1:Nmatrix
            for n = m:Nmatrix
                matrixout(m,n) = data(T,i);
                matrixout(n,m) = data(T,i);
                T = T+1;
            end
        end
    else
        for m = 1:Nmatrix-1
            for n = m+1:Nmatrix
                matrixout(n,m) = data(T,i);
                matrixout(m,n) = data(T,i);
                T = T+1;
            end
        end
    end
    lab_write_matrix(fullfile(Output_filepath,[Output_fileS '_' measures{i} '.txt']),matrixout);
end

end

function datainput = correctheader(datainput)
   if isnumeric(datainput)
       datainput = num2cell(datainput);
   end
   if isnumeric(datainput{1,1}) & isnumeric(datainput{2,1}) & isnumeric(datainput{1,2})
       for i = 1:size(datainput,2)
           header{1,i} = ['Measure' num2str(i)]; %#ok<AGROW>
       end
       for i = 1:size(datainput,1)
           measure{i,1} = ['Value' num2str(i)]; %#ok<AGROW>
       end
       measure = [{''};measure];
       datainput = [measure cat(1,header,datainput)];
   elseif isnumeric(datainput{2,1}) & ~isnumeric(datainput{1,2})
       for i = 1:size(datainput,1)-1
           measure{i,1} = ['Value' num2str(i)]; %#ok<AGROW>
       end
       measure = [{''};measure];
       datainput = [measure datainput];
   elseif isnumeric(datainput{1,2}) & ~isnumeric(datainput{2,1})
       for i = 1:size(datainput,2)-1
           header{1,i} = ['Measure' num2str(i)]; %#ok<AGROW>
       end
       header = [{''} header];
       datainput = cat(1,header,datainput);
   end
   tmp = find(cellfun(@isnumeric,datainput(1,2:end)));
   for i = tmp
       datainput{1,i+1} = num2str(datainput{1,i+1});
   end
   clearvars tmp
   tmp = find(cellfun(@isnumeric,datainput(2:end,1)));
   if isempty(tmp)
       return
   end
   for i = tmp'
       if isempty(datainput{i+1,1})
           datainput{i+1,1} = '';
       else
           datainput{i+1,1} = num2str(datainput{i+1,1});
       end
   end
   clearvars tmp
end