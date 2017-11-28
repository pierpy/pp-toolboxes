function lab_plot_distribution(data,skipsettings)
disp('Plot distribution')
if ~exist('skipsettings','var')
    skipsettings = false;
end

if ~exist('data','var') | ischar(data)
    if ~exist('data','var') | ~exist(data,'file')
        [data_file,data_filepath] = uigetfile('*.*','Select file');
    else
        [data_file,data_filepath] = lab_filename(data);
        clearvars data
    end
    data = lab_read_data(fullfile(data_filepath,data_file));
    cd(data_filepath);
end

if isnumeric(data) & skipsettings == false;
    settings.combine = true;
    settings = lab_set_plot_histogram(settings,data);
    if isempty(settings)
        return
    end
    
    if settings.combine == true
        data = data(:)';
    end
    if settings.excludezeros == true
        if size(data,1) > 1
            data(data == 0) = NaN;
        else
            data = setdiff(data,0);
            data = data(:);
        end
    end
elseif iscell(data)
    data = correctheader(data);
    Mtranspose = questdlg('Are subjects/trials in first row or column?','Subjects column/row','Cancel','Column','Row','Row');
    if strcmp(Mtranspose,'Column')
        data = data';
    end
    measures = data(2:end,1);
    data = cell2mat(data(2:end,2:end));
    
    settings.combine = false;
    settings = lab_set_plot_histogram(settings,data,measures);
    if isempty(settings)
        return
    end
    
    measures = measures(settings.selection,1);
    data = data(settings.selection,:);
    if settings.combine == true
        data = data(:)';
    end
    if settings.excludezeros == true
        if size(data,1) > 1
            data(data == 0) = NaN;
        else
            data = setdiff(data,0);
            data= data(:);
        end
    end
elseif isnumeric(data)
    data = data(:)';
else
    disp('Abort, wrong input data')
    return
end

f = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none');
if exist('data_filepath','var')
    [~,~,~,data_fileS] = lab_filename(data_file);
    data_fileS = regexprep(data_fileS,'_',' ');
    set(f,'Name',['Distribution ' data_fileS]);
else
    set(f,'Name','Distribution');
end
m1 = uimenu(f,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');

Nplots = size(data,1);
if ~exist('settings','var') | ~isfield(settings,'numvertical')
    if Nplots > 3
        NfigY = ceil(Nplots/4);
    else
        NfigY = 1;
    end
else
    NfigY = settings.numvertical;
end
NfigX = ceil(Nplots/NfigY);

for i = 1:Nplots
    subplot(NfigY,NfigX,i);
    datatmp = sort(data(i,:),'descend');
    plot(datatmp);
    if Nplots > 1 & exist('measures','var')
        title(regexprep(measures{i,1},'_',' '));
    end
    if exist('settings','var') & isfield(settings,'minval')
        set(gca,'YLim',[settings.minval settings.maxval]);
    end
end

end

function datainput = correctheader(datainput)
   if isnumeric(datainput{1,1}) & isnumeric(datainput{2,1}) & isnumeric(datainput{1,2})
       for i = 1:size(datainput,2)
           header{1,i} = ['Trial' num2str(i)];
       end
       for i = 1:size(datainput,1)
           measure{i,1} = ['Measure' num2str(i)];
       end
       measure = [{''};measure];
       datainput = [measure cat(1,header,datainput)];
   elseif isnumeric(datainput{2,1}) & ~isnumeric(datainput{1,2})
       for i = 1:size(datainput,1)-1
           measure{i,1} = ['Measure' num2str(i)];
       end
       measure = [{''};measure];
       datainput = [measure datainput];
   elseif isnumeric(datainput{1,2}) & ~isnumeric(datainput{2,1})
       for i = 1:size(datainput,2)-1
           header{1,i} = ['Trial' num2str(i)];
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