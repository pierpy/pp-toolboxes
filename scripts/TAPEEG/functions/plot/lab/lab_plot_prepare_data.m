function [DATA,PLOT,PlotList,Brain] = lab_plot_prepare_data(cfg,zerocolor)

if ~exist('zerocolor','var') | isempty(zerocolor)
    zerocolor = [1 1 1];
end

if ~isfield(cfg,'Store')
    cfg.Store = false;
end

DATA = [];
PLOT = [];
PlotList = [];
Brain = [];
if isfield(cfg,'Brain')
    Brain = cfg.Brain;
    % Exist if no valid brain information
    if ~isfield(Brain,'faces')
        disp('No valid brain template found')
        return
    end
    if ~isfield(Brain,'regions')
        Brain.regions = [];
    end
elseif isfield(cfg,'LOCS') & isempty(cfg.LOCS)
    disp('No valid LOCS information')
    return
end

% Set & correct data path
if isfield(cfg,'DATA_file') & cfg.DATA_file ~= 0 & ~isempty(cfg.DATA_file)
    if length(cfg.DATA_file) > 5
        [~,DATA_filepath,~,DATA_file] = lab_filename(cfg.DATA_file);
    end
    cfg.DATA_file = fullfile(DATA_filepath,regexprep(DATA_file,':','_'));
    cfg.DATA_file = fullfile(DATA_filepath,regexprep(DATA_file,';','_'));
    cfg.DATA_file = fullfile(DATA_filepath,regexprep(DATA_file,'//','_'));
    clearvars DATA_file DATA_filepath
elseif ~isempty('Brain')
    cfg.DATA_file = 'IS_Plot';
else
    cfg.DATA_file = 'Plot';
end

if isfield(cfg,'Mappings') & ~isempty(cfg.Mappings)
    LOCS = cfg.Mappings;
elseif isfield(cfg,'LOCS')
    LOCS = cfg.LOCS;
elseif isfield(cfg,'Brain')
    LOCS = cfg.Brain;
end
    
% skip processing if data not matching to brain template
cfg = lab_plot_check_numchans(cfg);
if cfg.Match == false
    DATA = [];
    PLOT = [];
    PlotList = [];
    return
end

% Sort data
[DATA,PLOT] = lab_plot_sort_data(cfg.DATA,cfg.PLOT);

% create color maps
for i = 1:length(PLOT)
    if isfield(PLOT(i),'ColorMode') & ~isempty(PLOT(i).ColorMode) & ~strcmp(PLOT(i).ColorMode,'color')
        PLOT(i).Color = lab_create_cmap(PLOT(i).ColorMode,zerocolor);
    elseif isfield(PLOT(i),'Color1') & ~isempty(PLOT(i).Color1)
        PLOT(i).Color = lab_create_cmap(PLOT(i).Color1,zerocolor);
    end
    if isfield(PLOT,'Threshold1') & isfield(PLOT,'Threshold2') & PLOT(i).Threshold1 > 0 & PLOT(i).Threshold2 > 0
        PLOT(i).Color = cat(1,lab_create_cmap(PLOT(i).Color2,zerocolor,32),lab_create_cmap(PLOT(i).Color1,zerocolor,32));
    end
end

% plot results
PlotList{1} = 1;
for i = 2:length(PLOT)
    if PLOT(i).AddPlot == true & (~isfield(PLOT(i),'labels') | PLOT(i).labels ~= 1)
        PlotList{i} = [PlotList{i-1} i]; %#ok<AGROW>
        PlotList{i-1} = []; %#ok<AGROW>
    else
        PlotList{i} = i; %#ok<AGROW>
    end
end
tmp = ~cellfun('isempty',PlotList);
PlotList = PlotList(tmp);
clearvars tmp

maxdata = length(DATA(:));
if length(PLOT) ~= maxdata
    DATA = [];
    return
end 
for datanr = 1:maxdata
    data = DATA(datanr).data;
    if isfield(PLOT(datanr),'Inverse') & PLOT(datanr).Inverse == true
        if PLOT(datanr).Threshold1 > 0
            T1 = data;
            T1(T1<=PLOT(datanr).Threshold1) = 0;
            T1(T1>0) = 1;
            if PLOT(datanr).Threshold2 > 0
                T2 = data;
                T2(T2<=PLOT(datanr).Threshold2) = 0;
                T2(T2>0) = 1;
            else
                T2 = ones(size(data));
            end
            data = zeros(size(data));
            if min(T2) == 0
                data(T2==0) = 0.5;
            end
            if min(T1) == 0
                data(T1==0) = 1;
            end
            clearvars T1 T2
            PLOT(datanr).MinValue = 0;
            PLOT(datanr).MaxValue = 1;
            Lidx = 1;
            if PLOT(datanr).Threshold2 > 0
                PLOT(datanr).Legend(Lidx,1).Color = PLOT(datanr).Color2;
                PLOT(datanr).Legend(Lidx,1).Text = [DATA(datanr).subject ' <= ' num2str(PLOT(datanr).Threshold2)];
                PLOT(datanr).Legend(Lidx,1).Text = regexprep(PLOT(datanr).Legend(Lidx,1).Text,'_',' ');
                PLOT(datanr).Legend(Lidx,1).Mode = 'Area';
                Lidx = Lidx + 1;
            end
            PLOT(datanr).Legend(Lidx,1).Color = PLOT(datanr).Color1;
            PLOT(datanr).Legend(Lidx,1).Text = [DATA(datanr).subject ' <= ' num2str(PLOT(datanr).Threshold1)];
            PLOT(datanr).Legend(Lidx,1).Text = regexprep(PLOT(datanr).Legend(Lidx,1).Text,'_',' ');
            PLOT(datanr).Legend(Lidx,1).Mode = 'Area';
            clearvars Lidx
        elseif min(data) >= 0 & max(data) <= 1
            data = 1 - data;
            PLOT(datanr).MinValue = 0;
            PLOT(datanr).MaxValue = 1;
            PLOT(datanr).Legend.Color = PLOT(datanr).Color;
            PLOT(datanr).Legend.Text = {num2str(PLOT(datanr).MinValue);num2str(PLOT(datanr).MaxValue)};
            PLOT(datanr).Legend.Mode = 'Colorbar';
        elseif min(data) < 0
            data = -data;
            tmp = PLOT(datanr).MaxValue;
            PLOT(datanr).MaxValue = -PLOT(datanr).MinValue;
            PLOT(datanr).MinValue = -tmp;
            clearvars tmp
            PLOT(datanr).Legend.Color = PLOT(datanr).Color;
            PLOT(datanr).Legend.Text = {num2str(PLOT(datanr).MinValue,3);num2str(PLOT(datanr).MaxValue,3)};
            PLOT(datanr).Legend.Mode = 'Colorbar';
        else
            data = 1/data;
            PLOT(datanr).MinValue = 0;
            PLOT(datanr).MaxValue = round(max(data));
            if PLOT(datanr).MaxValue < 1
                PLOT(datanr).MaxValue = 1;
            end
            PLOT(datanr).Legend.Color = PLOT(datanr).Color;
            PLOT(datanr).Legend.Text = {num2str(PLOT(datanr).MinValue);num2str(PLOT(datanr).MaxValue)};
            PLOT(datanr).Legend.Mode = 'Colorbar';
        end
    elseif isfield(PLOT,'Threshold1') & PLOT(datanr).Threshold1 > 0
        T1 = data;
        T1(T1>=PLOT(datanr).Threshold1) = 1;
        T1(T1<1) = 0;
        if PLOT(datanr).Threshold2 > 0
            T2 = data;
            T2(T2>=PLOT(datanr).Threshold2) = 1;
            T2(T2<1) = 0;
        else
            T2 = zeros(size(data));
        end
        data = zeros(size(data));
        if max(T2) == 1
            data(T2==1) = 0.5;
        end
        if max(T1) == 1
            data(T1==1) = 1;
        end
        PLOT(datanr).MinValue = 0;
        PLOT(datanr).MaxValue = 1;
        Lidx = 1;
        if PLOT(datanr).Threshold2 > 0
            PLOT(datanr).Legend(Lidx,1).Color = PLOT(datanr).Color2;
            PLOT(datanr).Legend(Lidx,1).Text = [DATA(datanr).subject ' >= ' num2str(PLOT(datanr).Threshold2)];
            PLOT(datanr).Legend(Lidx,1).Text = regexprep(PLOT(datanr).Legend(Lidx,1).Text,'_',' ');
            PLOT(datanr).Legend(Lidx,1).Mode = 'Area';
            Lidx = Lidx + 1;
        end
        PLOT(datanr).Legend(Lidx,1).Color = PLOT(datanr).Color1;
        PLOT(datanr).Legend(Lidx,1).Text = [DATA(datanr).subject ' >= ' num2str(PLOT(datanr).Threshold1)];
        PLOT(datanr).Legend(Lidx,1).Text = regexprep(PLOT(datanr).Legend(Lidx,1).Text,'_',' ');
        PLOT(datanr).Legend(Lidx,1).Mode = 'Area';
        clearvars Lidx
    elseif isfield(PLOT(datanr),'Channels') & ~isempty(PLOT(datanr).Channels)
        if isfield(PLOT(datanr),'MaxDistance') & PLOT(datanr).MaxDistance > 0 & isfield(LOCS,'x')
            distances = lab_distance(LOCS);
            distances = distances / min(distances(distances > 0));
            sigma = (PLOT(datanr).MaxDistance^2/(-2 * log(PLOT(datanr).WeightMaxDistance/100)))^0.5;
            for i = PLOT(datanr).Channels
                for j = 1:size(LOCS.x,2);
                    if distances(i,j) <= PLOT(datanr).MaxDistance
                        weights(i,j) = exp(-distances(i,j)^2 / (sigma^2 * 2)); %#ok<AGROW>
                    else
                        weights(i,j) = 0; %#ok<AGROW>
                    end
                end
            end
            data = max(weights,[],1);
            DATA(datanr).name = [DATA(datanr).name '_s' num2str(round(sigma*100)/100)];
            clearvars i j distances sigma weights
        else
            data = zeros(1,length(data));
            data(PLOT(datanr).Channels) = 1;
        end
        data = data*2 - 1;
        data(data==-1) = -1.2;
        PLOT(datanr).MinValue = -1.2;
        PLOT(datanr).MaxValue = 1.2;
        PLOT(datanr).Legend.Color = PLOT(datanr).Color1;
        if isfield(LOCS,'labels')
            tmp = sprintf('%s_',LOCS.labels{PLOT(datanr).Channels});
        elseif isfield(LOCS,'mappingstitle')
            tmp = sprintf('%s_',LOCS.mappingstitle{PLOT(datanr).Channels});
        else
            tmp = ['Channel' num2str(datanr) ' '];
        end
        PLOT(datanr).Legend.Text = tmp(1:end-1);
        PLOT(datanr).Legend.Text = regexprep(PLOT(datanr).Legend.Text,'_',' ');
        PLOT(datanr).Legend.Mode = 'Area';
        clearvars tmp
    else
        PLOT(datanr).Legend.Color = PLOT(datanr).Color;
        PLOT(datanr).Legend.Text = {num2str(PLOT(datanr).MinValue,3);num2str(PLOT(datanr).MaxValue,3)};
        PLOT(datanr).Legend.Mode = 'Colorbar';
    end
    DATA(datanr).data = data;
end

end