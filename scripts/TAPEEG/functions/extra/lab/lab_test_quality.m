% Test quality of spectral analysis using the file quality.xls
%
% written by F. Hatz Vumc 2013

function [Result,settings] = lab_test_quality

Result = [];

[data,header,outcome,~,settings] = lab_read_statistics([],1,0,0,0,1);
if isempty(data)
    return
else
    outcome = outcome';
    data = data';
end

if ~isfield(settings,'QUALITY')
    settings.QUALITY = [];
end
settings.QUALITY = lab_get_QUALITY(settings.QUALITY,1,0,0,0);

if size(outcome,1) > 1 & size(outcome,2) > 1
    disp('Abort, wrong input data');
    return
elseif size(outcome,1) > 1
    
end
if length(header.vars) ~= size(data,1)
    clearvars vars
end

Variables = {'cogfreq' 'Median Frequency';'peakfreq' 'Peak Frequency'; ...
    'peakamp' 'Peak Amplitude';'peak2min' 'Peak2Min';'peakratio' 'Peak Ratio'; ...
    'areapower' 'BandPower PeakArea';'bplast' 'BandPower High Frequency';[] ' '};
if exist('header','var') & isfield(header,'vars')
    Matching = header.vars(:);
else
    for i = 1:size(data,1)
        Matching{i,1} = ['Var' num2str(i,'%03d')];
    end
end
if size(data,1) <= 8
    Matching = [Matching Variables(1:size(data,1),2)];
else
    tmp = repmat(cellstr(''),size(data,1),1);
    tmp(1:8,1) = Variables(:,2);
    Matching = [Matching tmp];
end
Matching = lab_table_dialog(Matching,{'Name' 'Assigned'},'Matching Variables',0,{'char',Variables(:,2)'});
if isempty(Matching)
    return
end
pause(0.2);
Mtmp = Matching(:,2);
clearvars Matching
for i = 1:length(Mtmp)
    Matching{i} = Variables{strcmp(Variables(:,2),Mtmp{i}),1};
end
MatchingName = Mtmp(~strcmp(Mtmp,' '),1);

for i = 1:size(data,2)
    for j = 1:length(Matching)
        if ~isempty(Matching{j})
            R(i).(Matching{j}) = data(j,i);
        end
    end
end
clearvars data

settingstmp = create_settings(settings);

Result = calculate_result(R,settingstmp,outcome);
[~,Max]=max(Result.Kappa);
Result.MaxKappa = Result.Kappa(Max(1));
Result.Settings = settingstmp(Max);
Result.SettingsAll = settingstmp;
Result.Data.R = R;
Result.Data.Outcome = outcome;
Result.TestValues = settings.QUALITY;

data = [];
for i = 1:length(Result.SettingsAll)
    data = [data cell2mat(struct2cell(Result.SettingsAll(i)))];
end
data = data(1:end-1,:);
data = cat(1,data,Result.Kappa);

Vars = struct2cell(settings.QUALITY);
VarsNr = find(~cellfun(@isempty,(Vars(2:end,1))))+1;
Vars = Vars(VarsNr,1);
if size(data,1) == 4 & size(Vars,1) == 3
    for B = 1:2
        Rsum = zeros(length(Vars{1}),length(Vars{2}),length(Vars{3}));
        KappaTMP = reshape(Result.Kappa,[length(Vars{1}) length(Vars{2}) length(Vars{3})]);
        for i = 1+B:length(Vars{1})-B
            for j = 1+B:length(Vars{2})-B
                for m = 1+B:length(Vars{3})-B
                    tmp = KappaTMP(i-B:i+B,j-B:j+B,m-B:m+B);
                    Rsum(i,j,m) = mean(tmp(:));
                end
            end
        end
        if B == 1
            [~,Max3]=max(Rsum(:));
            Result.Settings3 = Result.SettingsAll(Max3);
        else
            [~,Max5]=max(Rsum(:));
            Result.Settings5 = Result.SettingsAll(Max5);
        end
        clearvars Rsumi j m tmp
    end
    clearvars B
elseif size(data,1) == 3 & size(Vars,1) == 2
    for B = 1:2
        Rsum = zeros(length(Vars{1}),length(Vars{2}));
        KappaTMP = reshape(Result.Kappa,[length(Vars{1}) length(Vars{2})]);
        for i = 1+B:length(Vars{1})-B
            for j = 1+B:length(Vars{2})-B
                tmp = KappaTMP(i-B:i+B,j-B:j+B);
                Rsum(i,j) = mean(tmp(:));
            end
        end
        if B == 1
            [~,Max3]=max(Rsum(:));
            Result.Settings3(Max3);
        else
            [~,Max5]=max(Rsum(:));
            Result.Settings5(Max5);
        end
        clearvars Rsumi j m tmp
    end
    clearvars B
end

if isfield(settings,'file')
    Filename = [settings.file 'Result.xlsx'];
    warning off
    mkdir(fullfile(settings.path,'TestQuality'));
    warning on
    Filepath = fullfile(settings.path,'TestQuality');
else
    [Filename,Filepath] = uiputfile('*.xlsx','Select File to store results');
    if Filename == 0
        return
    end
end

xlsout = Variables(VarsNr-1,2);
xlsout{end+1,1} = 'Kappa';
xlsout = [xlsout num2cell(data(:,Max))];
lab_write_xls(fullfile(Filepath,Filename),xlsout);

if exist('Max3','var')
    xlsout = Variables(VarsNr-1,2);
    xlsout{end+1,1} = 'Kappa';
    xlsout = [xlsout num2cell(data(:,Max3))];
    lab_write_xls(fullfile(Filepath,[Filename(1:end-5) '3.xlsx']),xlsout);
    
    xlsout = Variables(VarsNr-1,2);
    xlsout{end+1,1} = 'Kappa';
    xlsout = [xlsout num2cell(data(:,Max5))];
    lab_write_xls(fullfile(Filepath,[Filename(1:end-5) '5.xlsx']),xlsout);
end

xlsout = Variables(VarsNr-1,2);
xlsout{end+1,1} = 'Kappa';
xlsout = [xlsout num2cell(data)];
lab_write_xls(fullfile(Filepath,[Filename(1:end-5) '_All.xlsx']),xlsout');

xlsout = header.subjects';
xlsout = cat(1,xlsout,num2cell(Result.Matrix(Max,:)));
lab_write_xls(fullfile(Filepath,[Filename(1:end-5) '_Outcomes.xlsx']),xlsout');

fig = figure('Visible','off','Color',[1 1 1]);
plot(sort(Result.Kappa,'descend'));
title(['Kappa distribution (highest value: ' num2str(Result.Kappa(Max(1))) ')']);
lab_print_figure(fullfile(Filepath,[Filename(1:end-5) '.jpg']),fig);
close(fig);

if size(data,1) == 4 & size(Vars,1) == 3
    if settings.QUALITY.Npositiv == 2
        for B = 3:-1:1
            Varstmp = Vars(setdiff(1:3,B),1);
            VarsNames = Variables(VarsNr(setdiff(1:3,B))-1,2);
            settings2.QUALITY = settings.QUALITY;
            settings2.QUALITY.(Variables{VarsNr(B)-1}) = [];
            settingstmp = create_settings(settings2);
            tmp = calculate_result(R,settingstmp,outcome);
            tmp = reshape(tmp.Kappa,length(Varstmp{1}),length(Varstmp{2}));
            if length(Varstmp{1}) > 1 & length(Varstmp{2}) > 1
                fig = figure('Visible','off','Color',[1 1 1]);
                surface(Varstmp{1},Varstmp{2},tmp');
                colorbar
                title([VarsNames{1} ' vs. ' VarsNames{2}]);
                lab_print_figure(fullfile(Filepath,[Filename(1:end-5) '_' num2str(4-B) '.jpg']),fig);
                close(fig);
                eval(['Result.Kappa' num2str(4-B) '=tmp(:);']);
            end
        end
    elseif settings.QUALITY.Npositiv == 1
        for B = 1:3
            Varstmp = Vars(B,1);
            VarsNames = Variables(VarsNr(B)-1,2);
            settings2.QUALITY = settings.QUALITY;
            tmp = setdiff(1:3,B);
            settings2.QUALITY.(Variables{VarsNr(tmp(1))-1}) = [];
            settings2.QUALITY.(Variables{VarsNr(tmp(2))-1}) = [];
            settingstmp = create_settings(settings2);
            tmp = calculate_result(R,settingstmp,outcome);
            fig = figure('Visible','off','Color',[1 1 1]);
            plot(Vars{1},tmp.Kappa)
            title([VarsNames{1} ' (Kappa)']);
            lab_print_figure(fullfile(Filepath,[Filename(1:end-5) '_' num2str(B) '.jpg']),fig);
            close(fig);
            eval(['Result.Kappa' num2str(B) '=tmp.Kappa;']);
        end
    end
end

end

function settingstmp = create_settings(settings)

settingstmp.cogfreq = [];
settingstmp.peakfreq = [];
settingstmp.peakamp = [];
settingstmp.peak2min = [];
settingstmp.areapower = [];
settingstmp.bplast = [];
if ~isempty(settings.QUALITY.cogfreq)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.cogfreq)
        for j = 1:length(tmp)
            tmp(j).cogfreq = settings.QUALITY.cogfreq(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
if ~isempty(settings.QUALITY.peakfreq)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.peakfreq)
        for j = 1:length(tmp)
            tmp(j).peakfreq = settings.QUALITY.peakfreq(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
if ~isempty(settings.QUALITY.peakamp)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.peakamp)
        for j = 1:length(tmp)
            tmp(j).peakamp = settings.QUALITY.peakamp(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
if ~isempty(settings.QUALITY.peak2min)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.peak2min)
        for j = 1:length(tmp)
            tmp(j).peak2min = settings.QUALITY.peak2min(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
if ~isempty(settings.QUALITY.peakratio)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.peakratio)
        for j = 1:length(tmp)
            tmp(j).peakratio = settings.QUALITY.peakratio(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
if ~isempty(settings.QUALITY.areapower)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.areapower)
        for j = 1:length(tmp)
            tmp(j).areapower = settings.QUALITY.areapower(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
if ~isempty(settings.QUALITY.bplast)
    tmp = settingstmp;
    settingstmp = [];
    for i = 1:length(settings.QUALITY.bplast)
        for j = 1:length(tmp)
            tmp(j).bplast = settings.QUALITY.bplast(i);
        end
        settingstmp = [settingstmp tmp];
    end
end
for i = 1:length(settingstmp)
    settingstmp(i).Npositiv = settings.QUALITY.Npositiv;
end

end

function Result = calculate_result(R,settings,outcome)
   for i = 1:length(R)
       for j = 1:length(settings)
           Matrix(j,i) = lab_evaluate_quality2(R(i),settings(j));
       end
   end
   for i = 1:size(Matrix,1)
       classes = unique(outcome(:));
       classes2 = unique(Matrix(i,:)');
       datatmp = zeros(2,2);
       for j = 1:size(outcome,2)
           x = find(classes == outcome(1,j));
           y = find(classes2 == Matrix(i,j));
           if ~isempty(x) & ~isempty(y)
               datatmp(x,y) = datatmp(x,y) + 1;
           end
       end
       Result.KappaAll(i) = kappa(datatmp,0,0.05,1);
       Result.Kappa(i) = Result.KappaAll(i).kappa;
       Result.Correlation(i) = corr(Matrix(i,:)',outcome');
   end
   Result.Matrix = Matrix;
end