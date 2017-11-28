% Plot results of BrainWave, permutation or connectivity matrix
% in signal space
%
% lab_plot_elec(data,LOCS)
%
% 'DATA' and 'LOCS' are optional
%
% Written by F. Hatz Vumc Amsterdam 10/2012
% (distrubution only with permission of the author)

function settings = lab_plot_elec(DATA,LOCS)

global Main_Path

cfg = [];
if ~exist('DATA','var')
    DATA = [];
else
    if isstruct(DATA) & isfield(DATA,'x')
        cfg.LOCS = DATA;
        DATA = [];
    elseif isnumeric(DATA)
        DATA2 = DATA;
        clearvars DATA
        if size(DATA2,1) > size(DATA2,2)
            DATA2 = DATA2';
        end
        if size(DATA2,1) == size(DATA2,2)
            for i = 1:size(DATA2,3)
                D1(i,:) = diag(DATA2(:,:,i)); %#ok<AGROW>
            end
            D2 = lab_extract_tril_wodiag(DATA2)';
            DATA = [];
            for i = 1:size(D2,1)
                DATA(end+1,1).data = D2(i,:); %#ok<AGROW>
                DATA(end,1).name = ['Data' num2str(length(DATA))];
                DATA(end,1).measure = ['Data' num2str(length(DATA)) '_Connections'];
                DATA(end,1).subject = '';
                DATA(end,1).connections = true;
                DATA(end,1).nodes = false;
                DATA(end+1,1).data = D1(i,:); %#ok<AGROW>
                DATA(end,1).name = ['Data' num2str(length(DATA))];
                DATA(end,1).measure = ['Data' num2str(length(DATA)) '_Nodes'];
                DATA(end,1).subject = '';
                DATA(end,1).connections = false;
                DATA(end,1).nodes = true;
            end
        else
            for i = 1:size(DATA2,1)
                DATA(i,1).data = DATA2(1,:);
                DATA(i,1).name = ['Data' num2str(i)];
                DATA(i,1).measure = ['Data' num2str(i)];
                DATA(i,1).subject = '';
                DATA(end,1).connections = false;
                DATA(end,1).nodes = false;
            end
        end
        clearvars DATA2
    end
end
if exist('LOCS','var') & isstruct(LOCS) & isfield(LOCS,'x')
    cfg.LOCS = LOCS;
end

% Read electrodes file
if ~isfield(cfg,'LOCS')
    if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes.sfp'),'file')
        cfg.LOCS = fullfile(Main_Path,'electrodes.sfp');
    elseif ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes.els'),'file')
        cfg.LOCS = fullfile(Main_Path,'electrodes.els');
    elseif ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes.xyz'),'file')
        cfg.LOCS = fullfile(Main_Path,'electrodes.xyz');
    elseif exist(fullfile(pwd,'electrodes.sfp'),'file')
        cfg.LOCS = fullfile(pwd,'electrodes.sfp');
    elseif exist(fullfile(pwd,'electrodes.els'),'file')
        cfg.LOCS = fullfile(pwd,'electrodes.els');
    elseif exist(fullfile(pwd,'electrodes.xyz'),'file')
        cfg.LOCS = fullfile(pwd,'electrodes.xyz');
    else
        [cfg.LOCS,tmp] = uigetfile('*.els;*.sfp;*.xyz','Select electrodes file');
        if cfg.LOCS ~= 0
            cfg.LOCS = fullfile(tmp,cfg.LOCS);
        else
            return
        end
        clearvars tmp
    end
end
if ~isstruct(cfg.LOCS)
    if ischar(cfg.LOCS) & exist(cfg.LOCS,'file')
        cfg.loc_file = cfg.LOCS;
        cfg.LOCS = lab_read_locs(cfg.LOCS);
    else
        disp('Abort: no valid LOCS-file')
        return
    end
end

% set config and plot
cfg.DATA = DATA;
lab_set_plot_elec(cfg);

end
