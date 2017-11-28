% Collect results from microstates analysis
% Output is xls-files
%
% lab_collect_microstates(cfg)
%
% written by F. Hatz 2014

function lab_collect_microstates(cfg)
disp('Collect microstates data of all files')

skipprocessing = 0;
if ~exist('cfg','var')
    cfg = [];
    skipselection = false;
else
    skipselection = true;
end
Files = {};

if ~isfield(cfg,'CollectMicro') | ~isfield(cfg.CollectMicro,'searchfolder')
    [cfg,Files,skipprocessing] = lab_set_collect_microstates(cfg);
    if skipprocessing == 1
        return
    else
        pause(0.2);
    end
end

% turn log-file on
diary(fullfile(cfg.CollectMicro.searchfolder,'CollectMicrostates.log'));

% search files
if isempty(Files)
    Files = lab_collect_microstates_search(cfg,skipselection);
end
if isempty(Files) | skipprocessing == 1
    disp('no microstates results found')
    return
end

% Create Output-folder
warning off %#ok<WNOFF>
if ~isempty(cfg.CollectMicro.outputfolder)
    mkdir(fullfile(cfg.CollectMicro.searchfolder,cfg.CollectMicro.outputfolder));
    cfg.Output_filepath = fullfile(cfg.CollectMicro.searchfolder,cfg.CollectMicro.outputfolder);
else
    mkdir(fullfile(cfg.CollectMicro.searchfolder,'MicrostatesAnalysis'));
    cfg.Output_filepath = fullfile(cfg.CollectMicro.searchfolder,'MicrostatesAnalysis');
end
warning on %#ok<WNON>

Result = [];
for filenr = 1:size(Files,2)
    disp(['Collect Microstates data: ' lab_filename(Files{1,filenr})])
    R = load(Files{1,filenr});
    if isfield(R,'header') & isfield(R.header,'patient') & cfg.CollectMicro.subjectname == -1
        patient = R.header.patient;
    else
        [patient,cfg.CollectMicro,skipprocessing] = lab_subjectname(Files{1,filenr},cfg.CollectMicro);
    end
    if skipprocessing == 1
        return
    end
    if isfield(R,'Result') 
        if ~isfield(Result,'Nr')
            for i = 1:length(R.Result.Nr)
                Result(i).Nr = R.Result.Nr(i); %#ok<AGROW>
                Result(i).patient = cellstr(patient); %#ok<AGROW>
                Result(i).Template = R.Result.Template{i}; %#ok<AGROW>
                if isfield(R,'header')
                    Result(i).header = R.header; %#ok<AGROW>
                end
            end
        else
            for i = 1:length(Result)
                if length(R.Result.Template) >= i & size(R.Result.Template{i},2) == size(Result(i).Template,2) & ...
                        size(R.Result.Template{i},1) == size(Result(i).Template,1)
                    Result(i).patient = [Result(i).patient cellstr(patient)]; %#ok<AGROW>
                    Result(i).Template = cat(3,Result(i).Template,R.Result.Template{i}); %#ok<AGROW>
                end
            end
        end
    end
    if isfield(R,'Stats')
        for i = 1:length(R.Stats)
            if R.Stats(i).Nr == Result(i).Nr
                if ~isfield(Result,'Stats') | length(Result) < i | isempty(Result(i).Stats)
                    Result(i).Stats = {'';'Cluster';'Mean EV';'Mean Correlation'; ...
                        'Max Correlation';'Max Correlation Index';'Max GFP';'Max GFP Index'; ...
                        'Mean Duration';'Minimal Duration';'Maximal Duration'; ...
                        'Std Duration';'Sum Duration'}; %#ok<AGROW>
                    for j = 1:R.Stats(i).Nr
                        Result(i).Stats{end+1,1} = ['NextCluster' num2str(j)]; %#ok<AGROW>
                    end
                    for j = 1:R.Stats(i).Nr
                        Result(i).Stats{end+1,1} = ['PreviousCluster' num2str(j)]; %#ok<AGROW>
                    end
                end
                if ~isfield(Result,'StatsP') | isempty(Result(i).StatsP)
                    Result(i).StatsP = {'';'GEV';'KL';'MeanVariance';'CrossValidation'}; %#ok<AGROW>
                end
                tmp = {};
                for j = 1:R.Stats(i).Nr;
                    tmp{1,j} = [patient '_' num2str(j)]; %#ok<AGROW>
                end
                if ~isfield(R.Stats,'SumDuration')
                    R.Stats(i).SumDuration = zeros(1,R.Stats(i).Nr);
                end
                tmp = cat(1,tmp,num2cell(1:R.Stats(i).Nr),num2cell(R.Stats(i).MeanEV(:)'), ...
                    num2cell(R.Stats(i).MeanCorr(:)'),num2cell(R.Stats(i).MaxCorr(:)'), ...
                    num2cell(R.Stats(i).MaxCorrIdx(:)'),num2cell(R.Stats(i).MaxGfp(:)'), ...
                    num2cell(R.Stats(i).MaxGfpIdx(:)'),num2cell(R.Stats(i).MeanDuration(:)'), ...
                    num2cell(R.Stats(i).MinDuration(:)'),num2cell(R.Stats(i).MaxDuration(:)'), ...
                    num2cell(R.Stats(i).StdDuration(:)'),num2cell(R.Stats(i).SumDuration(:)'), ...
                    num2cell(R.Stats(i).PreviousCluster ./ repmat(sum(R.Stats(i).NextCluster,1),R.Stats(i).Nr,1)), ...
                    num2cell(R.Stats(i).NextCluster ./ repmat(sum(R.Stats(i).PreviousCluster,1),R.Stats(i).Nr,1)));
                Result(i).Stats = cat(2,Result(i).Stats,tmp); %#ok<AGROW>
                tmp = cat(1,patient,num2cell(R.Stats(i).GEV),num2cell(R.Stats(i).KL), ...
                    num2cell(R.Stats(i).MeanVariance),num2cell(R.Stats(i).CV));
                Result(i).StatsP = cat(2,Result(i).StatsP,tmp); %#ok<AGROW>
            end
        end
    end
end

% correct orientation
Result = lab_orient_microstates(Result);

settings = [];
StatsP = {};
for i = 1:length(Result)
    lab_write_xls(fullfile(cfg.Output_filepath,['Results_' num2str(Result(i).Nr) 'Clusters.xlsx']),Result(i).Stats);
    Output_filepath2 = fullfile(cfg.Output_filepath,[num2str(Result(i).Nr) 'Clusters']);
    warning off %#ok<WNOFF>
    mkdir(Output_filepath2)
    warning on %#ok<WNON>
    
    if isfield(Result(i),'header') & ~isempty(Result(i).header)
        headertopo = Result(i).header;
        headertopo.numtimeframes = Result(i).Nr;
        if isfield(cfg.CollectMicro,'dosef') & cfg.CollectMicro.dosef == true
            lab_write_sef(fullfile(Output_filepath2,[num2str(Result(i).Nr) 'Clusters.sef']),mean(Result(i).Template,3),headertopo);
        end
        if isfield(cfg.CollectMicro,'doplot') & cfg.CollectMicro.doplot == true & ...
            isfield(headertopo,'locs') & size(mean(Result(i).Template,3),1) == size(headertopo.locs.x,2) & ...
            (max(headertopo.locs.z) > 0 | min(headertopo.locs.z) < 0)
            disp(['   Write Microstate Topos (' num2str(Result(i).Nr) ' Clusters)'])
            if ~isfield(settings,'LOCS')
                settings.LOCS = headertopo.locs;
                settings.LOCS = lab_locs2sph(settings.LOCS);
            end
            PLOT.close = 0;
            PLOT.AddPlot = 0;
            for j = 1:size(mean(Result(i).Template,3),2)
                settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
                range = max(abs(mean(Result(i).Template(:,j,:),3)));
                PLOT.MinValue = -range;
                PLOT.MaxValue = range;
                settings = lab_plot_chans(mean(Result(i).Template(:,j,:),3),PLOT,settings);
                lab_print_figure(fullfile(Output_filepath2,['Cluster_' num2str(j) '.jpg']),settings.handleF);
                close(settings.handleF);
            end
            clearvars PLOT range f
        end
    end
    StatsP = cat(3,StatsP,Result(i).StatsP);
    StatsP{1,1,end} = ['Clust' num2str(Result(i).Nr)];
end

flag = true;
for i = 2:size(StatsP,2)
    tmp = unique(StatsP(1,i,:));
    if length(tmp) > 1
        flag = false;
    end
end
for i = 2:size(StatsP,1)
    tmp = unique(StatsP(i,1,:));
    if length(tmp) > 1
        flag = false;
    end
end
clearvars tmp
Result2 = [];
if flag == true
    tmp = StatsP(2:end,2:end,:);
    tmp = permute(tmp,[3 2 1]);
    TClust = StatsP(1,1,:);
    TClust = TClust(:);
    TVar = StatsP(2:end,1,1);
    for i = 1:size(tmp,3)
        xlsout = StatsP(1,:,1);
        xlsout = cat(1,xlsout,cat(2,TClust,tmp(:,:,i)));
        lab_write_xls(fullfile(cfg.Output_filepath,['Results_' TVar{i} '.xlsx']),xlsout);
        Result2.(TVar{i}) = xlsout;
    end
end

if exist('Result','var')
    save(fullfile(cfg.Output_filepath,'ResultsMicrostates.mat'),'Result','Result2');
end

disp('  Evaluate Microstates')
lab_evaluate_microstates(Result,Result2,cfg.Output_filepath);

% turn log-file off
diary off

end




