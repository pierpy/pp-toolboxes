% Collect bad channels from folder structure processed by TAPEEG
% Output is xls-file in main folder
%
% Result = lab_collect_badchannels(searchfolder)
%
% written by F. Hatz 2012

function Result = lab_collect_badchannels(searchfolder)

disp('Collect bad channels')

if exist('searchfolder','var') & exist(searchfolder,'dir')
    cfg.SEARCH.searchfolder = searchfolder;
end
cfg.SEARCH.searchstring = {'filt.info'};
[cfg,skipprocessing] = lab_set_searchstrings(cfg,1,0,{'strings','Verbose-file','Info-file','ICA-Container'});
if skipprocessing == 1;
    return
end
[FileList,cfg] = lab_search_files(cfg);
List = FileList.Filelist;
if isempty(List)
    disp('No files to process found')
    return
else
    List = List(:)';
end

Result.badchansSingleNames = {};
Result.badchansSingle = {};
Result.subject = {};
Subjectbadchans = [];
Subjectinterpolated = [];
Channels = {};
Numchannels = 0;
for N = 1:size(List,2)
    disp(['   read ' List{1,N}]);
    [~,~,Format] = lab_filename(List{1,N});
    if strcmp(Format,'vrb')
        read_vrb(List{1,N});
    elseif strcmp(Format,'mat')
        read_mat(List{1,N});
    else
        read_info(List{1,N});
    end
end
if ~isempty(Subjectbadchans)
    Result.Subjectbadchans = Subjectbadchans ./ repmat(Filesize,1,size(Subjectbadchans,2));
end
if ~isempty(Subjectinterpolated)
    Result.Subjectinterpolated = Subjectinterpolated ./ repmat(Filesize,1,size(Subjectinterpolated,2));
end

% Write results
if exist('Result','var')
    disp('   write results')
    Ftmp = pwd;
    warning off %#ok<WNOFF>
    mkdir(fullfile(cfg.SEARCH.searchfolder,'CollectBad'));
    warning on %#ok<WNON>
    cd(fullfile(cfg.SEARCH.searchfolder,'CollectBad'));
    
    resultfiles = [Result.filelist num2cell(Result.filelistbadchans)];
    resultfiles = cat(1,[cellstr('Files:') Channels'],resultfiles);
    if size(resultfiles,2) > 255
        filenameout = 'BadChannels_Files.xlsx';
    else
        filenameout = 'BadChannels_Files.xls';
    end
    lab_write_xls(filenameout,resultfiles);
    
    if isfield(Result,'filelistinterpolated')
        resultfilesI = [Result.filelist num2cell(Result.filelistinterpolated)];
        resultfilesI = cat(1,[cellstr('Files:') Channels'],resultfilesI);
        if size(resultfilesI,2) > 255
            filenameout = 'Interpolated_Files.xlsx';
        else
            filenameout = 'Interpolated_Files.xls';
        end
        lab_write_xls(filenameout,resultfilesI);
    end
    
    if isfield(Result,'Subjectbadchans')
        resultsubjects = [Result.subject num2cell(Result.Subjectbadchans)];
        resultsubjects = cat(1,[cellstr('Subjects:') Channels'],resultsubjects);
        if size(resultsubjects,2) > 255
            filenameout = 'BadChannels_Subjects.xlsx';
        else
            filenameout = 'BadChannels_Subjects.xls';
        end
        lab_write_xls(filenameout,resultsubjects);
    end
    
    if isfield(Result,'Subjectinterpolated')
        resultsubjectsI = [Result.subject num2cell(Result.Subjectinterpolated)];
        resultsubjectsI = cat(1,[cellstr('Subjects:') Channels'],resultsubjectsI);
        if size(resultsubjectsI,2) > 255
            filenameout = 'Interpolated_Subjects.xlsx';
        else
            filenameout = 'Interpolated_Subjects.xls';
        end
        lab_write_xls(filenameout,resultsubjectsI);
    end

    if isfield(Result,'badchansSingle')
        for i = 1:size(Result.badchansSingleNames,2)
            resultfiles = [Result.filelist(1:size(Result.badchansSingle{1,i},1)) num2cell(Result.badchansSingle{1,i})];
            resultfiles = cat(1,[cellstr('Files:') Channels'],resultfiles);
            if size(resultfiles,2) > 255
                filenameout = ['BadChannels_' Result.badchansSingleNames{1,i} '.xlsx'];
            else
                filenameout = ['BadChannels_' Result.badchansSingleNames{1,i} '.xls'];
            end
            filenameout = regexprep(filenameout,{':',',',';','\'},'');
            lab_write_xls(filenameout,resultfiles);
        end
    end
    
    cd(Ftmp);
end

    function read_vrb(Filename)
        [badchans,sizetmp,bad] = lab_read_badvrb(Filename);
        if min(badchans) ~= -1;
            [~,~,~,FilenameS] = lab_filename(Filename);
            [Subject,cfg] = lab_subjectname(Filename,cfg);
            Result.filelist{N,1} = FilenameS;
            if isempty(Result.subject) | max(strcmp(Result.subject,Subject)) == 0
                if Numchannels < max(badchans)
                    Numchannels = max(badchans);
                    Channels = cellstr(num2str((1:Numchannels)'));
                end
                Result.subject{end+1,1} = Subject;
                subjectnr = size(Result.subject,1);
                Filesize(subjectnr,1) = 0;
                Subjectbadchans(subjectnr,1:Numchannels) = zeros(1,Numchannels);
            else
                subjectnr = find(strcmp(Result.subject,Subject),1);
            end
            Filesize(subjectnr,1) = Filesize(subjectnr,1) + sizetmp;
            Result.filelistbadchans(N,1:Numchannels) = zeros(1,Numchannels);
            if badchans > 0
                Subjectbadchans(subjectnr,badchans) = Subjectbadchans(subjectnr,badchans) + sizetmp;
                Result.filelistbadchans(N,badchans) = 1;
            end
            for j = 1:size(bad.name,2)
                if strcmp(bad.name{1,j}(1:20),'Bad channels by file')
                    bad.name{1,j} = 'Bad channels by file';
                end
                Nname = find(strcmp(Result.badchansSingleNames,bad.name{1,j}),1);
                if isempty(Nname)
                    Result.badchansSingleNames = [Result.badchansSingleNames bad.name(1,j)];
                    Nname = size(Result.badchansSingleNames,2);
                    Result.badchansSingle{Nname} = zeros(N,Numchannels);
                end
                Result.badchansSingle{Nname}(N,1:Numchannels) = zeros(1,Numchannels);
                if isnan(bad.channels{1,j})
                    bad.channels{1,j} = [];
                end
                if ~isempty(bad.channels{1,j})
                    Result.badchansSingle{Nname}(N,bad.channels{1,j}) = 1;
                end
            end
        else
            disp(['   Skip ' lab_filename(Filename) ' - not valid'])
        end
        clearvars badchans
    end

    function read_mat(Filename)
        Mtmp = load(Filename);
        if isfield(Mtmp,'bad')
            badchans = Mtmp.bad;
        elseif isfield(Mtmp,'header') & isfield(Mtmp.header,'badchans')
            badchans = Mtmp.header.badchans;
        end
        if isfield(Mtmp,'header') & isfield(Mtmp.header,'interpolated')
            interpolated = Mtmp.header.interpolated;
        else
            interpolated = [];
        end
        if exist('badchans','var')
            if isfield(Mtmp.header,'numdatachannels')
                Numchannels2 = Mtmp.header.numdatachannels;
            else
                Numchannels2 = Mtmp.header.numchannels;
            end
            if Numchannels == 0
                Numchannels = Numchannels2;
                Channels = cellstr(Mtmp.header.channels(1:Numchannels,:));
            elseif Numchannels ~= Numchannels2
                disp(['   Skip ' lab_filename(Filename) ' - number of channels not matching'])
                return
            end
            [~,~,~,FilenameS] = lab_filename(Filename);
            [Subject,cfg] = lab_subjectname(Filename,cfg);
            Result.filelist{N,1} = FilenameS(1:end-4);
            if isempty(Result.subject) | max(strcmp(Result.subject,Subject)) == 0
                Result.subject{end+1,1} = Subject;
                subjectnr = size(Result.subject,1);
                Filesize(subjectnr,1) = 0;
                Subjectbadchans(subjectnr,1:Numchannels) = zeros(1,Numchannels);
                Subjectinterpolated(subjectnr,1:Numchannels) = zeros(1,Numchannels);
            else
                subjectnr = find(strcmp(Result.subject,Subject),1);
            end
            Filesize(subjectnr,1) = Filesize(subjectnr,1) + Mtmp.header.numtimeframes;
            Result.filelistbadchans(N,1:Numchannels) = zeros(1,Numchannels);
            if ~isempty(badchans) & badchans > 0
                Subjectbadchans(subjectnr,badchans) = Subjectbadchans(subjectnr,badchans) + Mtmp.header.numtimeframes;
                Result.filelistbadchans(N,badchans) = 1;
            end
            Result.filelistinterpolated(N,1:Numchannels) = zeros(1,Numchannels);
            if ~isempty(interpolated) & interpolated > 0
                Subjectinterpolated(subjectnr,interpolated) = Subjectinterpolated(subjectnr,interpolated) + Mtmp.header.numtimeframes;
                Result.filelistinterpolated(N,interpolated) = 1;
            end
        end
    end

    function read_info(Filename)
        [~,~,Format2,FilenameS] = lab_filename(Filename);
        if strcmp(Format2,'info')
            header = lab_read_eeginfo(Filename);
        else
            [~,header,cfg] = lab_read_data(Filename,cfg,1);
        end
        if isempty(header) | ~isfield(header,'badchans')
            disp(['   Skip ' lab_filename(Filename) ' - not valid'])
            return
        end
        if isfield(header,'numchannels')
            Numchannels2 = header.numchannels;
        elseif isfield(header,'badchans') & isfield(header,'goodchans') & isfield(header,'refchan') & isnumeric(header.ref_chan)
            Numchannels2 = max([header.badchans header.goodchans header.ref_chan]);
        elseif isfield(header,'badchans') & isfield(header,'goodchans')
            Numchannels2 = max([header.badchans header.goodchans]);
        end
        if Numchannels == 0
            Numchannels = Numchannels2;
            if isfield(header,'channels')
                Channels = cellstr(header.channels(1:Numchannels,:));
            else
                Channels = cellstr(num2str((1:Numchannels)'));
            end
        elseif Numchannels ~= Numchannels2
            disp(['   Skip ' lab_filename(Filename) ' - number of channels not matching'])
            return
        end
        [Subject,cfg] = lab_subjectname(Filename,cfg);
        Result.filelist{N,1} = FilenameS;
        if isempty(Result.subject) | max(strcmp(Result.subject,Subject)) == 0
            Result.subject{end+1,1} = Subject;
            subjectnr = size(Result.subject,1);
            Filesize(subjectnr,1) = 0;
            Subjectbadchans(subjectnr,1:Numchannels) = zeros(1,Numchannels);
            Subjectinterpolated(subjectnr,1:Numchannels) = zeros(1,Numchannels);
        else
            subjectnr = find(strcmp(Result.subject,Subject),1);
        end
        if ~isfield(header,'numtimeframes')
            header.numtimeframes = 1;
        end
        Filesize(subjectnr,1) = Filesize(subjectnr,1) + header.numtimeframes;
        Result.filelistbadchans(N,1:Numchannels) = zeros(1,Numchannels);
        if ~isempty(header.badchans) & header.badchans > 0
            Subjectbadchans(subjectnr,header.badchans) = Subjectbadchans(subjectnr,header.badchans) + header.numtimeframes;
            Result.filelistbadchans(N,header.badchans) = 1;
        end
        Result.filelistinterpolated(N,1:Numchannels) = zeros(1,Numchannels);
        if isfield(header,interpolated') & ~isempty(header.interpolated) & header.interpolated > 0
            Subjectinterpolated(subjectnr,header.interpolated) = Subjectinterpolated(subjectnr,header.interpolated) + header.numtimeframes;
            Result.filelistinterpolated(N,header.interpolated) = 1;
        end
    end
end
