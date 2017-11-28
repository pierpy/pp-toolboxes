% Function to xls-files with new montage info
%
% montage = lab_read_montage(cfg,MONT_file)
%                  or
% montage = lab_read_montage(MONT_file)
%
% written by F. Hatz 2012

function [montage,table] = lab_read_montage(MONT_file,cfg)


if ~exist('cfg','var')
    cfg = [];
end
if ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto')
    cfg.MAIN.auto = 0;
end
    
if ~exist('MONT_file','var')
    [MONT_file,MONT_filepath] = uigetfile('*.xls;*.xlsx','Select Montage-file');
    MONT_file = fullfile(MONT_filepath,MONT_file);
    clearvars MONT_filepath
end
if ~isempty(MONT_file)
    if ~exist('cfg','var') | ~isfield(cfg,'settings_path')
        cfg.settings_path = pwd;
    end
    if exist(fullfile(cfg.settings_path,MONT_file),'file')
        MONT_file = fullfile(cfg.settings_path,MONT_file);
    elseif exist(fullfile(cfg.settings_path,[MONT_file 'x']),'file')
        MONT_file = fullfile(cfg.settings_path,[MONT_file 'x']);
    elseif exist(fullfile(pwd,MONT_file),'file')
        MONT_file = fullfile(pwd,MONT_file);
    elseif exist(fullfile(pwd,[MONT_file 'x']),'file')
        MONT_file = fullfile(pwd,[MONT_file 'x']);
    end
end
if exist('MONT_file','var') & exist(MONT_file,'file')
    disp('   Read Montage from file')
    if ispc
        [~,~,mfile] = xlsread(MONT_file,1);
    else
        [~,~,mfile] = xlsread(MONT_file,1,'','basic');
    end
    clearvars MONT_file
end

if exist('mfile','var') & strcmp(mfile{1,1},'Montage')
    numchans = mfile{1,2};
    table = [];
    Mnr = 0;
    mfile = mfile(2:end,:);
elseif exist('mfile','var') & strcmp(mfile{1,2},'active') & strcmp(mfile{1,3},'reference')
    numchans = [];
    table = [];
    Mnr = 0;
end
if exist('Mnr','var')
    for i = 1:size(mfile,1)
        if strcmp(mfile{i,2},'active') & strcmp(mfile{i,3},'reference')
            Mnr = Mnr + 1;
            montage(1,Mnr).numchans = numchans; %#ok<AGROW>
            montage(1,Mnr).auxchans = 1; %#ok<AGROW>
            montage(1,Mnr).name = mfile{i,1}; %#ok<AGROW>
            montage(1,Mnr).chans = []; %#ok<AGROW>
            montage(1,Mnr).label = []; %#ok<AGROW>
        elseif Mnr > 0
            if ischar(mfile{i,2})
                if length(mfile{i,2}) > 3 & strcmp(mfile{i,2}(1:3),'AUX')
                    tmp{1,2} = 1;
                    tmp2 = str2num(mfile{i,2}(4:end)); %#ok<ST2NM>
                    if ~isempty(tmp2)
                        tmp{1,1} = tmp2;
                    else
                        tmp{1,1} = 0;
                        tmp{1,2} = 0;
                    end
                    clearvars tmp2
                else
                    tmp{1,2} = 0;
                    tmp2 = str2num(mfile{i,1}); %#ok<ST2NM>
                    if ~isempty(tmp2)
                        tmp{1,1} = tmp2;
                    else
                        tmp{1,1} = 0;
                    end
                    clearvars tmp2
                end
            else
                tmp{1,2} = 0;
                if ~isempty(mfile{i,2}) & isnumeric(mfile{i,2})
                    tmp{1,1} = mfile{i,2};
                else
                    tmp{1,1} = 0;
                end
            end
            if ischar(mfile{i,3})
                if length(mfile{i,3}) > 3 & strcmp(mfile{i,3}(1:3),'AUX')
                    tmp{1,4} = 1;
                    tmp2 = str2num(mfile{i,3}(4:end)); %#ok<ST2NM>
                    if ~isempty(tmp2)
                        tmp{1,3} = tmp2;
                    else
                        tmp{1,3} = 0;
                        tmp{1,4} = 0;
                    end
                    clearvars tmp2
                elseif length(mfile{i,3}) == 3 & strcmp(mfile{i,3},'AVG')
                    tmp{1,3} = mfile{i,3};
                    tmp{1,4} = 0;
                elseif length(mfile{i,3}) == 4 & strcmp(mfile{i,3},'LAPL')
                    tmp{1,3} = mfile{i,3};
                    tmp{1,4} = 0;
                else
                    tmp{1,4} = 0;
                    tmp2 = str2num(mfile{i,3}); %#ok<ST2NM>
                    if ~isempty(tmp2)
                        tmp{1,3} = tmp2;
                    else
                        tmp{1,3} = 0;
                    end
                    clearvars tmp2
                end
            else
                tmp{1,4} = 0;
                if ~isempty(mfile{i,3}) & isnumeric(mfile{i,3}) & ~isnan(mfile{i,3})
                    tmp{1,3} = mfile{i,3};
                else
                    tmp{1,3} = 0;
                end
            end
            if max(tmp{1,1}) > montage(1,Mnr).numchans
                montage(1,Mnr).numchans = max(tmp{1,1}); %#ok<AGROW>
            end
            if max(tmp{1,3}) > montage(1,Mnr).numchans
                montage(1,Mnr).numchans = max(tmp{1,3}); %#ok<AGROW>
            end
            if tmp{1,2} == 1 & max(tmp{1,1}) > montage(1,Mnr).auxchans
                montage(1,Mnr).auxchans = max(tmp{1,1}); %#ok<AGROW>
            end
            if tmp{1,4} == 1 & max(tmp{1,3}) > montage(1,Mnr).auxchans
                montage(1,Mnr).auxchans = max(tmp{1,3}); %#ok<AGROW>
            end
            montage(1,Mnr).chans = cat(1,montage(1,Mnr).chans,tmp); %#ok<AGROW>
            if isnumeric(mfile{i,1})
                tmp = cellstr(num2str(mfile{i,1}));
            else
                tmp = mfile(i,1);
            end
            montage(1,Mnr).label = cat(1,montage(1,Mnr).label,tmp); %#ok<AGROW>
            clearvars tmp
        end
    end
    if Mnr > 0
        table = lab_montage2table(montage);
    end
else
    disp('    Wrong file format (montage)')
    montage = [];
    table = [];
end