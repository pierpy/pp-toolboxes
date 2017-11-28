function lab_evaluate_microstates(Result,Result2,Filepath)

if ~exist('Result','var') | ~exist('Result2','var')
    [Filename,Filepath] = uigetfile({'.mat'},'Select Microstates-Results');
    if isempty(Filename) | isnumeric(Filename)
        return
    end
    load(fullfile(Filepath,Filename));
end
if ~exist('Result','var') | ~exist('Result2','var')
    disp('No result file of ''Collect Microstates'' selected');
    return
end

datKL = cell2mat(Result2.KL(2:end,2:end));
for i = 1:size(datKL,2);
    datKL(:,i) = (datKL(:,i) - min(datKL(:,i))) / (max(datKL(:,i)) - min(datKL(:,i)));
end
datCV = cell2mat(Result2.CrossValidation(2:end,2:end));
for i = 1:size(datCV,2);
    datCV(:,i) = (datCV(:,i) - min(datCV(:,i))) / (max(datCV(:,i)) - min(datCV(:,i)));
end
warning off %#ok<WNOFF>
for i = 1:size(datCV,2)
    % tmp = find(diff(datCV(:,i)) == max(diff(datCV(:,i))),1,'first');
    tmp = find(diff(datCV(:,i)) > 0,1,'first');
    if ~isempty(tmp)
        Opt(1,i) = tmp;
    else
        Opt(1,i) = size(datCV,1);
    end
    tmp = corner(1:size(datKL,1),datKL(:,i));
    if ~isempty(tmp)
        Opt(2,i) = tmp;
    else
        Opt(2,i) = size(datKL,1);
    end
end
warning on %#ok<WNON>

fig1 = figure('Color',[1 1 1],'NumberTitle','off','Name','Clustering -- CV','Visible','off');
Clusters = Result2.KL(2:end,1);
Subjects = Result2.KL(1,2:end);

plot(datCV);
title('Clustering -- Cross Validation')
set(gca,'XTick',1:length(Clusters),'XTickLabel',Clusters);
rotateXLabelsImage(gca,90);
lab_print_figure(fullfile(Filepath,'Microstates-CrossValidation.tif'),fig1);

plot(datKL);
title('Clustering -- Krzanowski-Lai')
set(gca,'XTick',1:length(Clusters),'XTickLabel',Clusters);
rotateXLabelsImage(gca,90);
lab_print_figure(fullfile(Filepath,'Microstates-Krzanowski-Lai.tif'),fig1);

dataout = [];
events.POS = [];
events.DUR= [];
events.OFF = [];
events.TYP = {};
for i = 1:size(Opt,2)
    events.POS = [events.POS int64(size(dataout,2)+1)];
    events.DUR = [events.DUR int64(size(Result(Opt(1,i)).Template(:,:,i),2))];
    events.OFF = [events.OFF int64(0)];
    events.TYP = [events.TYP Subjects(1,i)];
    dataout = cat(2,dataout,Result(Opt(1,i)).Template(:,:,i));
end
header = Result(1).header;
header = lab_reduce_header(header,1:size(dataout,1));
header.goodchans = 1:size(dataout,1);
header.badchans = [];
if isfield(header,'interpolated')
    header.interpolated = [];
end
header.events = events;
if isfield(header,'bad')
    header = rmfield(header,'bad');
end
lab_write_sef(fullfile(Filepath,'All_OptCV_Microstates'),dataout,header);
