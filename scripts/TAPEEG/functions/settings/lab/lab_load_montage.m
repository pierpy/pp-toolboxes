function Montage = lab_load_montage(Montage,cfg,header,Filename)

global Main_Path
    
if ~exist('Filename','var')
    Filename = 'Montage.xls';
end
if ~exist('header','var')
    header = [];
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('Montage','var')
    Montage = [];
end

if isempty(Montage)
    if isfield(cfg,'settings_path') & exist(fullfile(cfg.settings_path,Filename),'file')
        montage_file = fullfile(cfg.settings_path,Filename);
    elseif isfield(cfg,'settings_path') & exist(fullfile(cfg.settings_path,[Filename 'x']),'file')
        montage_file = fullfile(cfg.settings_path,[Filename 'x']);
    elseif exist(fullfile(pwd,Filename),'file')
        montage_file = fullfile(pwd,Filename);
    elseif exist(fullfile(pwd,[Filename 'x']),'file')
        montage_file = fullfile(pwd,[Filename 'x']);
    end
    if exist('montage_file','var')
        Montage = read_montage(montage_file,cfg,header);
        if ~isempty(Montage)
            Montage = lab_reduce_montage(Montage,cfg,header,true);
        end
    end
end

settings.montage = Montage;

List_montage = {};
if ~isempty(Main_Path) & exist(fullfile(Main_Path,'montage'),'dir')
    List_montage = lab_search(fullfile(Main_Path,'montage'),{'*.xls','*.xlsx'},true,true,1);
end
List_montage = cat(1,{'Select File'},List_montage(:));

if exist('montage_file','var')
    settings.montage_file = montage_file;
    List_montage = cat(1,{montage_file},List_montage(:));
else
    settings.montage_file = '';
    List_montage = cat(1,{''},List_montage(:));
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Montage-File','montage_file'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = List_montage;
Formats(end,1).callback =  {@read_montage,{'montage','montage_file'},'montage_file',cfg,header};
Formats(end,1).size = 270;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Montage','montage'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = [130 70];
Formats(end,1).enable = 'inactive';
Formats(end,1).span = [2 1];

Prompt(end+1,:) = {'Edit',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@edit_montage,'montage','montage',header};

Prompt(end+1,:) = {'Reduce',''};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@lab_reduce_montage,'montage','montage',cfg,header};

[settings,Cancelled] = inputsdlg(Prompt,'Load Montage',Formats,settings);
if Cancelled == 1
    Montage = [];
    return
end
Montage = settings.montage;

    function [montage,montage_file] = read_montage(montage_file,cfg,header)
        if ~isempty(montage_file) & strcmp(montage_file,'Select File')
            [montage_file,montage_path] = uigetfile({'*.xls;*.xlsx','Excel-File'},'Select Montage-File');
            montage_file = fullfile(montage_path,montage_file);
        end
        if isempty(montage_file) | ~exist(montage_file,'file')
            montage_file = '';
            montage = [];
            return
        end
        montage = lab_read_montage(montage_file);
        if isempty(montage)
            return
        end
        montage = lab_reduce_montage(montage,cfg,header,true);
    end

    function montage = edit_montage(montage,header)
        titles={'name','active','reference'};
        if ~isempty(montage) & isfield(montage,'chans')
            montagetmp = montage(1);
            montage = montagetmp;
            clearvars montagetmp
            table(:,1) = montage.label;
            for i = 1:size(table,1)
                if ischar(montage.chans{i,1})
                    table{i,2} = montage.chans{i,1};
                elseif montage.chans{i,2} == 1
                    table{i,2} = ['AUX' num2str(montage.chans{i,1})];
                else
                    table{i,2} = num2str(montage.chans{i,1});
                end
                if ischar(montage.chans{i,3})
                    table{i,3} = montage.chans{i,3};
                elseif montage.chans{i,4} == 1
                    table{i,3} = ['AUX' num2str(montage.chans{i,3})];
                else
                    table{i,3} = num2str(montage.chans{i,3});
                end
            end
        elseif isfield(header,'numdatachannels') & header.numdatachannels == 257
            montage.name = 'EDFexport';
            montage.numchans = 257;
            table={'AF3','AF4','AF7','AF8','C1','C2','C3','C4','C5','C6','CP1','Cp2', ...
                'CP3','CP4','CP5','CP6','CPz','Cz','F1','F10','F2','F3','F4','F5','F6', ...
                'F7','F8','F9','FC1','FC2','FC3','FC4','FC5','FC6','FCz','FP1','FP2', ...
                'FPz','FT10','FT7','FT8','FT9','Fz','O1','O2','Oz','P1','P10','P2','P3', ...
                'P4','P5','P6','P9','PO3','PO4','PO7','PO8','POz','Pz','T10','T3','T4', ...
                'T5','T6','T9','TP10','Tp7','TP8','TP9'; ...
                '34','12','46','10','44','185','59','183','64','194','79','143','66', ...
                '164','76','172','81','257','29','226','5','36','224','48','222', ...
                '47','2','252','24','207','42','206','49','213','15','37','18','26', ...
                '219','62','211','67','21','116','150','126','88','169','142','87','153', ...
                '86','162','106','109','140','97','161','119','101','210','69','202','96', ...
                '170','68','190','84','179','94'; ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190','94 190','94 190', ...
                '94 190','94 190','94 190','94 190','94 190','94 190'}';
        else
            montage.name = 'Montage';
            if isfield(header,'numdatachannels')
                montage.numchans = header.numdatachannels;
            else
                montage.numchans = [];
            end
            table={'1','1',''};
        end
        tinfo = lab_table_dialog(table,titles,'EDF channels',1);
        if ~isempty(tinfo)
            montagetmp = montage;
            montage = [];
            montage.name = montagetmp.name;
            montage.numchans = montagetmp.numchans;
            try
                for i = 1:size(tinfo,1)
                    if length(tinfo{i,2}) > 3 & strcmp(tinfo{i,2}(1:3),'AUX')
                        montage.chans{i,1} = str2num(tinfo{i,2}(4:end)); %#ok<ST2NM>
                        montage.chans{i,2} = 1;
                    else
                        tmp = str2num(tinfo{i,2}); %#ok<ST2NM>
                        if isempty(tmp) | isnan(tmp)
                            montage.chans{i,1} = tinfo{i,2};
                        else
                            montage.chans{i,1} = tmp;
                        end
                        clearvars tmp 
                        montage.chans{i,2} = 0;
                    end
                    if length(tinfo{i,3}) > 3 & strcmp(tinfo{i,3}(1:3),'AUX')
                        montage.chans{i,3} = str2num(tinfo{i,3}(4:end)); %#ok<ST2NM>
                        montage.chans{i,4} = 1;
                    else
                        tmp = str2num(tinfo{i,3}); %#ok<ST2NM>
                        if isempty(tmp) | isnan(tmp)
                            montage.chans{i,3} = tinfo{i,3};
                        else
                            montage.chans{i,3} = tmp;
                        end
                        clearvars tmp 
                        montage.chans{i,4} = 0;
                    end
                end
                montage.label = tinfo(:,1);
                if ~isfield(montage,'numchans') | isempty(montage.numchans)
                    maxchan = 0;
                    for i = 1:length(montage.chans)
                        if max(montage.chans{i}) > maxchan
                            maxchan = max(montage.chans{i});
                        end
                    end
                    montage.numchans = maxchan;
                end
            catch %#ok<CTCH>
                montage = [];
            end
        else
            montage = [];
        end
        clearvars titles table tinfo
    end
end