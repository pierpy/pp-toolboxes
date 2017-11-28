% Collect distances from IS-results
% Output is xls-files
%
% Result = lab_collect_distance(cfg)
%
% written by F. Hatz 2013

function Result = lab_collect_distances(cfg)

disp('Collect distance data of all files')
Result = [];
Files = {};

if ~exist('cfg','var')
    cfg = [];
    skipselection = false;
else
    skipselection = true;
end

if ~isfield(cfg,'CollectDistance') | isempty(cfg.CollectDistance)
    [cfg,Files,skipprocessing] = lab_set_collect_distances(cfg);
    if skipprocessing == 1
        return
    end
end

% search files
if isempty(Files)
    Files = lab_collect_distances_search(cfg,skipselection);
end

if isempty(Files)
    disp('No Distance-Files found')
    return
end

if ~isfield(cfg.CollectDistance,'Threshold') | isempty(cfg.CollectDistance.Threshold)
    cfg.CollectDistance.Threshold = 5;
end

if isempty(cfg.CollectGraph.outputfolder)
    disp('Abort: specify output folder')
    return
end
warning off %#ok<WNOFF>
mkdir(fullfile(cfg.CollectDistance.searchfolder,cfg.CollectDistance.outputfolder));
Output_path = fullfile(cfg.CollectDistance.searchfolder,cfg.CollectDistance.outputfolder);
warning on %#ok<WNON>

% collect data
Resulttmp = [];
Patients = [];
for filenr = 1:length(Files)
    MAT = load(Files{filenr});
    if isfield(MAT,'result')
        if filenr == 1
            for i = 1:size(MAT.result.file,2)
                tmp = strfind(MAT.result.file{i},'_');
                if ~isempty(tmp)
                    list{i} = MAT.result.file{i}(tmp(end)+1:end); %#ok<AGROW>
                else
                    list{i} = MAT.result.file{i}; %#ok<AGROW>
                end
                tmp = strfind(list{i},'.');
                if ~isempty(tmp)
                    list{i} =list{i}(1:tmp(end)-1); %#ok<AGROW>
                end
            end
        end
        Resulttmp = cat(4,Resulttmp,MAT.result.distance);
        if isfield(MAT,'patient') & cfg.CollectDistance.subjectname == -1
            Patients = [Patients cellstr(MAT.patient)]; %#ok<AGROW>
        elseif ~isempty(cfg.CollectDistance.subjectname)
            Patients = [Patients cellstr(lab_subjectname(Files{filenr},cfg.CollectDistance))]; %#ok<AGROW>
        else
            Patients = [Patients cellstr(lab_subjectname(Files{filenr}))]; %#ok<AGROW>
        end
    end
    clearvars MAT
end

if exist('Resulttmp','var') & ~isempty(Resulttmp)
    save(fullfile(Output_path,'Resultstmp.mat'),'Resulttmp','patients','list');
    Resulttmp = permute(Resulttmp,[1 2 4 3]);
    for i = 1:length(list)
        calc = Resulttmp(:,:,:,i);
        eval(['Result.' list{i} '=calc;']);
        calc(Resulttmp(:,:,:,i) > cfg.CollectDistance.Threshold) = 0;
        calc(Resulttmp(:,:,:,i) <= cfg.CollectDistance.Threshold) = 1;
        if size(calc,2) > 255
            xlsfileout = fullfile(Output_path,['Distance_' list{i} '.xlsx']);
        else
            xlsfileout = fullfile(Output_path,['Distance_' list{i} '.xls']);
        end
        lab_write_xls(xlsfileout,mean(calc,3));
        f = figure('Visible','off');
        colormap('gray');
        cmap = colormap;
        colormap(flipud(cmap));
        imagesc(mean(calc,3));
        lab_print_figure(fullfile(Output_path,['Distance_' list{i} '.tif']),f);
    end
    if ~isempty(Patients)
        Result.patients = Patients;
    end
    save(fullfile(Output_path,'Results.mat'),'Result');
end

return