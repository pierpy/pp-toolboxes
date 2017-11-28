% Calculate CoV
%
% [Result,settings] = lab_CoV(data,settings)
%
% data               = array (subjects x measures)
% settings.bootstrap = number of samples drawn (0 = off)
%
% Written by F. Hatz 2014

function [Result,settings] = lab_CoV(data,settings)

if ~exist('settings','var')
    settings = [];
end
if ~exist('data','var') | isempty(data)
    [data,header,~,~,settings] = lab_read_statistics([],-1,0,1,1,1);
    if isempty(data)
        Result = [];
        return
    end
end
if ~exist('header','var')
    [File,Path] = uiputfile('*.xls','Select File to store results');
    header.file = File;
    header.path = Path;
end

if isempty(settings) | ~isfield(settings,'bootstrap')
    [settings,skipprocessing] = lab_set_CoV(settings);
    if skipprocessing == 1
        Result = [];
        return
    end
end

% Create Output-Folder
warning off; %#ok<WNOFF>
mkdir(fullfile(header.path,'CoV'));
warning on; %#ok<WNON>
settings.path = fullfile(header.path,'CoV');

% Correct NaN-values
exclude = find(sum(isnan(data),2)>0);
if ~isempty(exclude)
    if ~isfield(settings,'button')
        settings.button = questdlg('NaN-values: Exclude or Replace by zero','NaN-values','Exclude','Replace','Replace');
    end
    if strcmp(settings.button,'Exclude')
        disp(['Exclude Subjects/Trials with NaN-values (' num2str(exclude) ')'])
        include = 1:size(data,1);
        include = setdiff(include,exclude);
        data = data(include,:,:);
    else
        disp('Replace NaN-values by 0')
        data(isnan(data)) = 0;
    end
end

if isfield(header,'vars')
    Result.Measure = header.vars;
else
    for i = 1:size(data,2)
        Result.Measure{1,i} = ['Rating' num2str(i)];
    end
    clearvars i
end

% CoV
Result.CoV = std(data,[],1) ./ mean(data,1);

% ICC jackknife SE
if isfield(settings,'ci') & strcmp(settings.ci,'jackknife')
    disp('Calculate jackknife SE')
    ICC_JK = lab_jackknife(data,@do_CoV);
    Result.jackknife_low = ICC_JK.value_low;
    Result.jackknife_high = ICC_JK.value_high;
end

% CoV bootstrapping
if isfield(settings,'ci') & strcmp(settings.ci,'bootstrap') & settings.bootstrap > 0
    fprintf('Calculate bootstrapping')
    rng('default');
    rng('shuffle');
    settings.counter = 0;
    settings.dcounter = floor(settings.bootstrap / 20);
    [CoV_ci,CoV_B] = bootci(settings.bootstrap,@do_CoV,data);
    Result.CoV_Bootstrapping = CoV_B;
    Result.CoV_Bootstrapping_CI = CoV_ci;
    disp(':')
end

% Write Results
if isfield(header,'file')
    xlsout = {'';'Coefficient of variation'};
    xlsout(1,2:length(Result.Measure)+1) = Result.Measure;
    xlsout(2,2:end) = num2cell(Result.CoV);
    if isfield(Result,'jackknife_low')
        xlsout(3:4,1) = {'Jackknife SE low';'Jackknife SE high'};
        xlsout(3,2:end) = num2cell(Result.jackknife_low);
        xlsout(4,2:end) = num2cell(Result.jackknife_high);
    end
    if isfield(Result,'CoV_Bootstrapping')
        xlsout2(1:11,1) = {'';['Bootstrapping (' num2str(settings.bootstrap) 'x)']; ...
            'Mean';'Std';'Median';'Min';'Max';'95% CI-Low';'95% CI-High'; ...
            'Quartile-Low';'Quartile-High'};
        for i = 1:size(CoV_B,2)
            xlsout2{3,i+1} = mean(CoV_B(:,i));
            xlsout2{4,i+1} = std(CoV_B(:,i));
            xlsout2{5,i+1} = median(CoV_B(:,i));
            xlsout2{6,i+1} = min(CoV_B(:,i));
            xlsout2{7,i+1} = max(CoV_B(:,i));
            xlsout2{8,i+1} = CoV_ci(1,i);
            xlsout2{9,i+1} = CoV_ci(2,i);
            tmp = sort(CoV_B(:,i));
            xlsout2{10,i+1} = tmp(ceil(0.25*length(tmp)));
            xlsout2{11,i+1} = tmp(ceil(0.75*length(tmp)));
        end
        xlsout = cat(1,xlsout,xlsout2);
    end
    lab_write_xls(fullfile(settings.path,[header.file '_CoV.xlsx']),xlsout);
end

% Plot Bootstrapping
if isfield(Result,'CoV_Bootstrapping')
    lab_write_xls(fullfile(settings.path,[header.file '_CoV_bootstrap.csv']),Result.CoV_Bootstrapping);
        
    fig1 = figure('NumberTitle','off','Name','CoV-Bootstrap','MenuBar','none','Color',[1 1 1]);
    m1 = uimenu(fig1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig1,'Label','Edit');
    uimenu(m2,'Label','15 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,15));
    uimenu(m2,'Label','30 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,30));
    uimenu(m2,'Label','45 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,45));
    uimenu(m2,'Label','60 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,60));
    uimenu(m2,'Label','90 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,90));
    boxplot(CoV_B,'labels',Result.Measure,'Colors','k','Symbol','+k');
    title('CoV (Bootstrap)');
    lab_print_figure(fullfile(settings.path,[header.file '_CoV_bootstrap.tif']),fig1);
end

    function CoV = do_CoV(data)
        CoV = std(data,[],1) ./ mean(data,1);
        if isfield(settings,'counter')
            settings.counter = settings.counter + 1;
            if mod(settings.counter,settings.dcounter) == 0
                fprintf('.');
            end
        end
    end

end

function set_X_Angle(H1,H2,labels,angle)
   text_h = findobj(gca, 'Type', 'text');
   for cnt = 1:length(text_h)
       tmp = labels{length(labels)-cnt+1};
       if length(tmp) > 5
           set(text_h(cnt),'Rotation',angle, ...
               'String',tmp, ...
               'HorizontalAlignment', 'right');
       end
   end
end
