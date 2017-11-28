% Collect and sort files in folderstructure
%
% lab_collect_filedata(cfg)
%
% written by F. Hatz 2013

function lab_collect_filedata(cfg)

if ~exist('cfg','var')
    cfg = [];
end

[cfg,skipprocessing] = lab_set_collect_filedata(cfg);
if skipprocessing == 1
    return
end

List = lab_search(cfg.CollectFiles.searchfolder,cfg.CollectFiles.stringfilename);
if isempty(List)
    disp('no files found')
    return
end
List = List(:);
for i = 1:size(List,1)
    [List{i,1},List{i,2}] = lab_filename(List{i,1});
end
if size(cfg.CollectFiles.stringpath,1) > 0
    good = [];
    for i = 1:size(cfg.CollectFiles.stringpath,1)
        tmp = ~cellfun('isempty',strfind(List(:,2),cfg.CollectFiles.stringpath{i,1}));
        tmp = find(tmp);
        good = union(good,tmp);
        clearvars tmp
    end
    if isempty(good)
        disp('no files found')
        return
    else
        List = List(good,:);
    end
end

good = 1:size(List,1);
for i = 1:size(cfg.CollectFiles.estringfilename,1)
    tmp = ~cellfun('isempty',strfind(List(:,1),cfg.CollectFiles.estringfilename{i,1}));
    tmp = find(tmp);
    good = setdiff(good,tmp);
    clearvars tmp
end
if isempty(good)
    disp('no files found')
    return
end

List = List(good,:);
good = 1:size(List,1);
for i = 1:size(cfg.CollectFiles.estringpath,1)
    tmp = ~cellfun('isempty',strfind(List(:,2),cfg.CollectFiles.estringpath{i,1}));
    tmp = find(tmp);
    good = setdiff(good,tmp);
    clearvars tmp
end
if isempty(good)
    disp('no files found')
    return
end

List = List(good,:);
clearvars good
if size(List,2) > 0
    disp ('Select Files')
    for i = 1:size(List,1)
        strlist{i,1} = fullfile(List{i,2},List{i,1}); %#ok<AGROW>
    end
    strdefault = 1:size(List,1);
    selection = listdlg('PromptString','Files:','SelectionMode','multiple', ...
        'ListString',strlist,'InitialValue',strdefault,'CancelString','None','ListSize',[450 400]);
    if ~isempty(selection)
        List = List(selection,:);
    else
        return
    end
    clearvars selection strlist strdefault
    pause(0.2);
end

if ~isfield(cfg.CollectFiles,'outputfolder') | isempty(cfg.CollectFiles.outputfolder);
    cfg.CollectFiles.outputfolder = fullfile(cfg.CollectFiles.searchfolder,'CollectFiles');
end
if ~exist(cfg.CollectFiles.outputfolder,'dir')
    warning off %#ok<WNOFF>
    mkdir(cfg.CollectFiles.outputfolder);
    warning on %#ok<WNON>
end

if cfg.CollectFiles.docopy == true
    cd(cfg.CollectFiles.outputfolder);
    for i = 1:size(List,1)
        copyfile(fullfile(List{i,2},List{i,1}),cfg.CollectFiles.outputfolder,'f');
    end
end

if cfg.CollectFiles.doaverage == true
    if isfield(cfg.CollectFiles,'grandaverage') & cfg.CollectFiles.grandaverage == true
        [~,~,~,FilenameS] = lab_filename(List{1,1});
        cfg.CollectFiles.averagesearchstrings{1,1} = FilenameS;
        if isempty(cfg.CollectFiles.averagestring)
            cfg.CollectFiles.averagestring = {'GrandAverage'};
        end
        for i = 1:size(List,1)
            ListAVG{1,i} = fullfile(List{i,2},List{i,1}); %#ok<AGROW>
        end
    else
        for i = 1:size(cfg.CollectFiles.averagesearchstrings,1)
            tmp = ~cellfun('isempty',strfind(List(:,1),cfg.CollectFiles.averagesearchstrings{i,1}));
            ListAVGtmp = List(tmp,:);
            clearvars tmp
            for j = 1:size(ListAVGtmp,1)
                tmp = strfind(ListAVGtmp{j,1},cfg.CollectFiles.averagesearchstrings{i,1});
                AVGsearch = [ListAVGtmp{j,1}(1:tmp-1) ListAVGtmp{j,1}(tmp+length(cfg.CollectFiles.averagesearchstrings{i,1}):end)];
                clearvars tmp
                if i == 1
                    ListAVG{j,i} = fullfile(ListAVGtmp{j,2},ListAVGtmp{j,1}); %#ok<AGROW>
                    ListAVGsearch{j,1} = AVGsearch; %#ok<AGROW>
                else
                    tmp = ~cellfun('isempty',strfind(ListAVGsearch(:,1),AVGsearch));
                    if length(tmp) == 1
                        ListAVG{tmp,i} = fullfile(ListAVGtmp{j,2},ListAVGtmp{j,1}); %#ok<AGROW>
                    end
                    clearvars tmp
                end
            end
        end
    end
    if exist('ListAVG','var')
        for i = 1:size(ListAVG,1)
            data = [];
            for j = 1:size(ListAVG,2)
                if ~isempty(ListAVG{i,j})
                    [datatmp,header] = lab_read_data(ListAVG{i,j});
                    if ~isempty(cfg.CollectFiles.INTERPOLATE)
                        [datatmp,header,cfg] = do_interpolate_bad(datatmp,header,cfg);
                    end
                    if isempty(data)
                        data = datatmp;
                    elseif size(data,1) == size(datatmp,1)
                        tmp = min(size(data,2),size(datatmp,2));
                        data = cat(3,data(:,1:tmp,:),datatmp(:,1:tmp));
                        header.numtimeframes = tmp;
                        clearvars tmp
                    end
                end
            end
            if size(data,3) > 1
                data = mean(data,3);
            end
            outputfile = lab_filename(regexprep(ListAVG{i,1},cfg.CollectFiles.averagesearchstrings{1,1},cfg.CollectFiles.averagestring{1,1}));
            lab_save_data(data,header,fullfile(cfg.CollectFiles.outputfolder,outputfile));
            ListAVG2{i,1} = outputfile; %#ok<AGROW>
            List2{i,1} = outputfile; %#ok<AGROW>
            List2{i,2} = cfg.CollectFiles.outputfolder; %#ok<AGROW>
            clearvars outputfile
        end
        ListAVG = cat(2,ListAVG2,ListAVG);
        if size(ListAVG,2) > 255
            fileout = fullfile(cfg.CollectFiles.outputfolder,'AverageFiles.xlsx');
        else
            fileout = fullfile(cfg.CollectFiles.outputfolder,'AverageFiles.xls');
        end
        lab_write_xls(fileout,ListAVG);
        List = List2;
    else
        disp('Abort: no files to average')
        return
    end
end

if cfg.CollectFiles.dogfp == true
    xlsout = {'Subject','GFP'};
    for i = 1:size(List,1)
        [datatmp,header,cfg] = lab_read_data(fullfile(List{i,2},List{i,1}),cfg);
        if ~isempty(datatmp)
            if ~isempty(cfg.CollectFiles.INTERPOLATE)
                [datatmp,header,cfg] = do_interpolate_bad(datatmp,header,cfg);
            end
            [patient,cfg.CollectFiles,skipprocessing] = lab_subjectname(fullfile(List{i,2},List{i,1}),cfg.CollectFiles);
            if skipprocessing == 1
                return
            end
            tmp = find(strcmp(xlsout(:,1),patient));
            if isempty(tmp)
                xlsout{end+1,1} = patient; %#ok<AGROW>
                tmp = size(xlsout,1);
            end
            [~,~,gfp] = lab_compute_gfp(datatmp,header);
            if ~isempty(gfp)
                xlsout{tmp,2} = [xlsout{tmp,2} gfp];
            end
        end
    end
    for i = 2:size(xlsout,1)
        if ~isempty(xlsout{i,2})
            xlsout{i,2} = mean(xlsout{i,2});
        else
            xlsout{i,2} = 0;
        end
    end
    disp('Write GFP.xls')
    fileout = fullfile(cfg.CollectFiles.outputfolder,'GFP.xls');
    lab_write_xls(fileout,xlsout);
end

end

function [datatmp,header,cfg] = do_interpolate_bad(datatmp,header,cfg)
   if ~isfield(header,'badchans')
       header.badchans = [];
   end
   if isfield(cfg.CollectFiles.INTERPOLATE,'definebad') & ~isempty(cfg.CollectFiles.INTERPOLATE.definebad)
       header.badchans = cfg.CollectFiles.INTERPOLATE.definebad;
   end
   if isfield(cfg.CollectFiles.INTERPOLATE,'BAD') & ~isempty(cfg.CollectFiles.INTERPOLATE.BAD)
       badchans = lab_detect_bad(datatmp,header,cfg.CollectFiles.INTERPOLATE.BAD);
       header.badchans = union(header.badchans,badchans);
       if size(header.badchans,1) > 1
           header.badchans = header.badchans';
       end
   end
   if isfield(cfg.CollectFiles.INTERPOLATE,'dointerpolate') & cfg.CollectFiles.INTERPOLATE.dointerpolate == true
       [datatmp,header] = lab_interpolate_bad(datatmp,header);
   end
end
