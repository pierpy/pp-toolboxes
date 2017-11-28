% Generate subject name from file-/folder-name:
% By defining how many underscores in filename are included or (with minus
% sign) wich level of foldername, the subject name is defined.
%
% [subjectname,cfg,skipprocessing] = lab_subjectname(Filename,cfg)
%
% written by F. Hatz 2013

function [subjectname,cfg,skipprocessing] = lab_subjectname(Filename,cfg)

skipprocessing = 0;
[Output,Foldername,Filename,Mainfolder] = lab_prepare_subjectname(Filename);
tmp = strfind(Filename,'_');
if length(tmp) > 1
    SecondVar = Filename(tmp(1)+1:tmp(2)-1);
elseif length(tmp) == 1
    SecondVar = Filename(tmp(1)+1:end);
else
    SecondVar = '';
end
clearvars tmp

if ~exist('cfg','var')
    if length(SecondVar) > 4 & strcmp(SecondVar(1:5),'Micro')
        cfg.subjectname = 1;
    else
        cfg.subjectname = 0;
    end
elseif ~isfield(cfg,'subjectname')
    if ~isempty(Output{1}) & isfield(cfg,'settings_path') & length(Mainfolder) > 1 & ...
            (strcmp(cfg.settings_path,Mainfolder) | strcmp(cfg.settings_path,Mainfolder(1:end-1)))
        cfg.subjectname = -1;
    elseif length(SecondVar) > 4 & strcmp(SecondVar(1:5),'Micro')
        cfg.subjectname = 1;
    else
        cfg.subjectname = 0;
    end
    if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
        Prompt = {};
        Formats = {};
        Prompt(end+1,:) = {'Number of underscores in subject name',''};
        Formats(end+1,1).type = 'text';
        if ~isempty(Output{1})
            Prompt(end+1,:) = {[Output{1} ' ' Output{2}],'subjectname'};
        else
            Prompt(end+1,:) = {Output{2},'subjectname'};
        end
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).limits = [-99 99];
        Formats(end,1).size = 40;
        
        [cfg,Cancelled] = inputsdlg(Prompt,'Subject name',Formats,cfg);
        if Cancelled == 1
            skipprocessing = 1;
            if ~isempty(Output{1})
                subjectname = get_name(Foldername,1,1);
            else
                subjectname = Filename;
            end
            return
        else
            pause(0.2);
        end
    end
end

% Create subject name
if cfg.subjectname < 0
    subjectname = get_name(Foldername,-cfg.subjectname,1);
else
    subjectname = get_name(Filename,cfg.subjectname+1,0);
end
cfg.patient = subjectname;

end


function NameOut = get_name(NameIn,Index,Inverse)
    tmp = union(strfind(NameIn,'_'),strfind(NameIn,' '));
    if size(tmp,1) > 1
        tmp = tmp';
    end
    if Inverse == 0
        tmp = [tmp,length(NameIn)+1];
        if Index <= length(tmp)
            if tmp(Index) > 1
                NameOut = NameIn(1:tmp(Index)-1);
            else
                NameOut = NameIn;
            end
        else
            NameOut = NameIn;
        end
    else
        tmp = [0 tmp];
        if length(tmp) == 1
            NameOut = NameIn;
        elseif Index > (length(tmp))
            NameOut = NameIn(1:tmp(2)-1);
        else
            tmp2 = length(tmp)+1-Index;
            if tmp2 >= 1 & tmp2 < length(tmp)
                NameOut = NameIn(tmp(tmp2)+1:tmp(tmp2+1)-1);
            elseif tmp2 >= 1 & tmp2 == length(tmp)
                NameOut = NameIn(tmp(tmp2)+1:end);
            else
                NameOut = NameIn(1:tmp(2)-1);
            end
        end
    end
end


