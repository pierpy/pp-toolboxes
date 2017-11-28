% Calculate ICC
%
% [Result,settings] = lab_ICC(data,settings)
%
% data               = array (measures x ratings)
% settings.mode      = '1-1' '1-k' A-1' 'A-k' 'C-1' 'C-k' (mode for ICC)
% settings.bootstrap = number of samples drawn (0 = off)
% settings.alpha     = alpha-value for ICC calculation
%
% ICC - Code from Mathworks Fileexchange
%
% Written by F. Hatz 2013

function [Result,settings] = lab_ICC(data,settings)

if ~exist('settings','var')
    settings = [];
end
if ~exist('data','var') | isempty(data)
    settings.collectfiles = 1;
    [data,header,~,~,settings] = lab_read_statistics(settings,-1,0,1,1,1);
    if isempty(data)
        Result = [];
        return
    end
end
if ~exist('header','var')
    [File,Path] = uigetfile('*.xls','Select File to store results');
    header.file = File;
    header.path = Path;
end

if isempty(settings) | ~isfield(settings,'cse')
    [settings,skipprocessing] = lab_set_ICC(settings);
    if skipprocessing == 1
        Result = [];
        return
    end
end

if ~isfield(settings,'clustervars') | settings.clustervars == 1
    Prompt = {'Number of Judgements per measure','clustervars'};
    Formats.type = 'edit';
    Formats.format = 'integer';
    Formats.limits = [2 9999999];
    Formats.size = 30;
    [settings,Cancelled] = inputsdlg(Prompt,'Number of Judgements per measure',Formats,settings);
    if isempty(settings) | Cancelled == 1 | settings.clustervars <= 1
        Result = [];
        return
    else
        pause(0.2);
        settings.numclusters = floor(size(data,2) / settings.clustervars);
    end
end
if ~isfield(settings,'numclusters')
    settings.numclusters = floor(size(data,2) / settings.clustervars);
end

% Create Output-Folder
warning off; %#ok<WNOFF>
mkdir(fullfile(header.path,'ICC'));
warning on; %#ok<WNON>
settings.path = fullfile(header.path,'ICC');

% Reshape data
data = reshape(data,[size(data,1),settings.clustervars,settings.numclusters]);
if ~strcmp(settings.orientation,'measures')
    data = permute(data,[3 2 1]);
end

% Correct for NaN
exclude = [];
for i = 1:size(data,3)
    exclude = union(exclude,find(sum(isnan(data(:,:,i)),2)>0));
end
if ~isempty(exclude)
    if ~isfield(settings,'button')
        settings.button = questdlg('NaN-values: Exclude or Replace by zero','NaN-values','Exclude','Replace','Replace');
    end
    if strcmp(settings.button,'Exclude')
        if strcmp(settings.orientation,'measures')
            disp(['Exclude Subjects/Trials with NaN-values (' num2str(exclude(:)') ')'])
        else
            disp(['Exclude Measures with NaN-values (' num2str(exclude(:)') ')'])
        end
        include = 1:size(data,1);
        include = setdiff(include,exclude);
        data = data(include,:,:);
    else
        disp('Replace NaN-values by 0')
        data(isnan(data)) = 0;
    end
end

% Create names for measures and variables
if isfield(header,'vars')
    if strcmp(settings.orientation,'measures')
        varstmp = reshape(header.vars,[size(data,2),size(data,3)]);
    
        vars = {};
        for i = 1:size(varstmp,1)
            tmp = strfind(varstmp{i,1},'_');
            if ~isempty(tmp)
                vars{1,i} = regexprep(varstmp{i,1}(tmp(end)+1:end),{'*','/','_'},' '); %#ok<AGROW>
            else
                vars{1,i} = regexprep(varstmp{i,1},{'*','/','_'},' '); %#ok<AGROW>
            end
        end
        Result.Variables = vars;
        clearvars vars
        
        vars = {};
        for i = 1:size(varstmp,2)
            tmp = strfind(varstmp{1,i},'_');
            if ~isempty(tmp)
                vars{1,i} = regexprep(varstmp{1,i}(1:tmp(end)-1),{'*','/','_'},' '); %#ok<AGROW>
            else
                vars{1,i} = regexprep(varstmp{1,i},{'*','/','_'},' '); %#ok<AGROW>
            end
        end
        Result.Measure = vars;
        clearvars vars
    else
        varstmp = reshape(header.vars,[size(data,2),size(data,1)]);
        vars = {};
        for i = 1:size(varstmp,1)
            tmp = strfind(varstmp{i,1},'_');
            if ~isempty(tmp)
                vars{1,i} = regexprep(varstmp{i,1}(tmp(end)+1:end),{'*','/','_'},' '); %#ok<AGROW>
            else
                vars{1,i} = regexprep(varstmp{i,1},{'*','/','_'},' '); %#ok<AGROW>
            end
        end
        Result.Variables = vars;
        clearvars vars
        
        Result.Measure = header.subjects(:)';
    end
    clearvars varstmp i
else
    for i = 1:size(data,3)
        Result.Measure{1,i} = ['Measure' num2str(i)];
    end
    for i = 1:size(data,2)
        Result.Variables{1,i} = ['Rating' num2str(i)];
    end
    clearvars i
end

% Plot input data
if settings.plotdata == true
    disp('Plot input data')
    for Nvar = 1:size(data,3)
        fig1 = figure('NumberTitle','off','Name','ICC','MenuBar','none','Color',[1 1 1],'Visible','off');
        m1 = uimenu(fig1,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
        plot(data(:,:,Nvar),'.');
        if length(Result.Variables) < 5
            legend(Result.Variables,'Location','NorthOutside','Orientation','horizontal');
        else
            legend(Result.Variables,'Location','NorthOutside');
        end
        xlabel(['ICC (' settings.mode '): ' Result.Measure{1,Nvar}]);
        if isfield(header,'file')
            lab_print_figure(fullfile(settings.path,[settings.file  '_' Result.Measure{1,Nvar} '_ICCdata.tif']),fig1);
            close(fig1);
        else
            set(fig1,'Visible','on');
        end
    end
end

settings.clustervars = size(data,2);
settings.numclusters = size(data,3);
data = reshape(data,size(data,1),size(data,2) * size(data,3));

% Calculate ICC
disp('Calculate ICC')
[ICC_R,ICC_LB,ICC_UB] = do_ICC(data);
Result.ICC = ICC_R;
Result.LowerBorder = ICC_LB;
Result.UpperBorder = ICC_UB;

% ICC jackknife SE
if isfield(settings,'ci') & strcmp(settings.ci,'jackknife')
    disp('Calculate jackknife SE')
    ICC_JK = lab_jackknife(data,@do_ICC);
    Result.jackknife_low = ICC_JK.value_low;
    Result.jackknife_high = ICC_JK.value_high;
end

% ICC bootstrapping
if isfield(settings,'ci') & strcmp(settings.ci,'bootstrap') & settings.bootstrap > 0
    fprintf('Calculate bootstrapping')
    rng('default');
    rng('shuffle');
    settings.counter = 0;
    settings.dcounter = floor(settings.bootstrap / 100);
    [ICC_ci,ICC_B] = bootci(settings.bootstrap,@do_ICC,data);
    Result.ICC_Bootstrapping = ICC_B;
    Result.ICC_Bootstrapping_CI = ICC_ci;
    disp(':')
end

% Write Results
if isfield(header,'file')
    xlsout = {settings.mode;'ICC'; ...
        ['Lower border (alpha:' num2str(settings.alpha) ')']; ...
        ['Upper border (alpha:' num2str(settings.alpha) ')']};
    xlsout(1,2:length(Result.Measure)+1) = Result.Measure;
    xlsout(2,2:end) = num2cell(Result.ICC);
    xlsout(3,2:end) = num2cell(Result.LowerBorder);
    xlsout(4,2:end) = num2cell(Result.UpperBorder);
    if isfield(Result,'jackknife_low')
        xlsout(5:6,1) = {'Jackknife SE low';'Jackknife SE high'};
        xlsout(5,2:end) = num2cell(Result.jackknife_low);
        xlsout(6,2:end) = num2cell(Result.jackknife_high);
    end
    if isfield(Result,'ICC_Bootstrapping')
        xlsout2(1:10,1) = {['Bootstrapping (' num2str(settings.bootstrap) 'x)']; ...
            'ICC-mean';'Std';'ICC-Median';'Min';'Max';'95% CI-Low';'95% CI-High'; ...
            'Quartile-Low';'Quartile-High'};
        for i = 1:size(ICC_B,2)
            xlsout2{2,i+1} = mean(ICC_B(:,i));
            xlsout2{3,i+1} = std(ICC_B(:,i));
            xlsout2{4,i+1} = median(ICC_B(:,i));
            xlsout2{5,i+1} = min(ICC_B(:,i));
            xlsout2{6,i+1} = max(ICC_B(:,i));
            xlsout2{7,i+1} = ICC_ci(1,i);
            xlsout2{8,i+1} = ICC_ci(2,i);
            tmp = sort(ICC_B(:,i));
            xlsout2{9,i+1} = tmp(ceil(0.25*length(tmp)));
            xlsout2{10,i+1} = tmp(ceil(0.75*length(tmp)));
        end
        xlsout = cat(1,xlsout,xlsout2);
    end
    lab_write_xls(fullfile(settings.path,[header.file '_ICC.xlsx']),xlsout');
end

% Plot Bootstrapping
if isfield(Result,'ICC_Bootstrapping')
    lab_write_xls(fullfile(settings.path,[header.file '_ICC_bootstrap.csv']),Result.ICC_Bootstrapping);

    fig2 = figure('NumberTitle','off','Name','ICC-Bootstrap','MenuBar','none','Color',[1 1 1]);
    m1 = uimenu(fig2,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig2,'Label','Edit');
    uimenu(m2,'Label','15 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,15));
    uimenu(m2,'Label','30 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,30));
    uimenu(m2,'Label','45 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,45));
    uimenu(m2,'Label','60 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,60));
    uimenu(m2,'Label','90 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,90));
    boxplot(Result.ICC_Bootstrapping,'labels',Result.Measure,'Colors','k','Symbol','+k');
    title(['ICC (' settings.mode ' - Bootstrap)']);
    set_X_Angle([],[],Result.Measure,45);
    if isfield(header,'file')
        lab_print_figure(fullfile(settings.path,[header.file '_ICC_bootstrap.tif']),fig2);
    end
elseif isfield(Result,'jackknife_low') & isfield(Result,'jackknife_high')
    dataout = cat(1,Result.jackknife_low,Result.jackknife_low,Result.ICC,Result.jackknife_high,Result.jackknife_high);
    fig2 = figure('NumberTitle','off','Name','ICC-Jackknife','MenuBar','none','Color',[1 1 1]);
    m1 = uimenu(fig2,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(fig2,'Label','Edit');
    uimenu(m2,'Label','15 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,15));
    uimenu(m2,'Label','30 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,30));
    uimenu(m2,'Label','45 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,45));
    uimenu(m2,'Label','60 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,60));
    uimenu(m2,'Label','90 grad','Callback',@(~,~)set_X_Angle([],[],Result.Measure,90));
    boxplot(dataout,'labels',Result.Measure,'Colors','k','Symbol','+k');
    title('CronbachAlpha - Jackknife');
    set_X_Angle([],[],Result.Measure,45);
    if isfield(header,'file')
        lab_print_figure(fullfile(settings.path,[header.file '_ICC_jackknife.tif']),fig2);
    end
end

    function [ICC_R,ICC_LB,ICC_UB] = do_ICC(data)
        datatmp = reshape(data,[size(data,1),settings.clustervars,settings.numclusters]);
        for j = 1:size(datatmp,3)
            [ICC_R(1,j),ICC_LB(1,j),ICC_UB(1,j)] = ICC(datatmp(:,:,j),settings.mode,settings.alpha); %#ok<AGROW>
        end
        Idx = find(isnan(ICC_R));
        if ~isempty(Idx)
            ICC_R(Idx) = nanmean(ICC_R);
            ICC_LB(Idx) = nanmean(ICC_LB);
            ICC_UB(Idx) = nanmean(ICC_UB);
        end
        if isfield(settings,'counter')
            settings.counter = settings.counter + 1;
            if mod(settings.counter,settings.dcounter) == 0
                fprintf('.');
            end
        end
    end

end


function set_X_Angle(H1,H2,labels,angle)
   text_h = findobj(gca,'Type','text');
   for cnt = 1:length(text_h)
       tmp = labels{length(labels)-cnt+1};
       if length(tmp) > 4
           set(text_h(cnt),'Rotation',angle, ...
               'String',tmp,'HorizontalAlignment','right');
       end
   end
end
