function [lag,cfg] = lab_get_egi_lag(samplingrate,Filepath,cfg)
    
global Main_Path
if ~exist('cfg','var')
    cfg = [];
end

if exist(fullfile(Filepath,'info1.xml'),'file')
    Info1 = xmltools(fullfile(Filepath,'info1.xml'));
    HWfilter = false;
    if isfield(Info1,'child')
        for i = 1:size(Info1.child,2)
            if isfield(Info1.child(1,i),'tag') & ...
                    strcmp(Info1.child(1,i).tag,'HARDWAREFILTERADJUSTED') & ...
                    strcmp(Info1.child(1,i).value,'true')
                HWfilter = true;
            end
            if isfield(Info1.child(1,i),'child')
                for j = 1:size(Info1.child(1,i).child,2)
                    if isfield(Info1.child(1,i).child(1,j),'tag') & ...
                            strcmp(Info1.child(1,i).child(1,j).tag,'HARDWAREFILTERADJUSTED') & ...
                            strcmp(Info1.child(1,i).child(1,j).value,'true')
                        HWfilter = true;
                    end
                end
            end
        end    
    end
    if HWfilter == true
        disp('    lag of events already corrected while recording')
        lag = 0;
        return
    end
end


InfoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info,'info.xml',Filepath);
amplifier = char(InfoObj.getAmpSerialNumber);
if isempty(amplifier)
    amplifier = 'DefaultEmpty';
end

LAGS = [];
if isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'LAGS')
    LAGS = cfg.EXTRA.LAGS;
elseif exist(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'),'file')
    load(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'))
elseif exist(fullfile(Main_Path,'MFF_Lags.xls'),'file')
    tmp = lab_read_xls(fullfile(Main_Path,'MFF_Lags.xls'),1);
    if ~isempty(tmp) & strcmp(tmp{1,1},'Frequency')
        for i = 2:size(tmp,1)
            for j = 2:size(tmp,2)
                if ~isempty(tmp{i,j}) & ~isnan(tmp{i,j})
                    LAGS.(tmp{i,1}).(['S' num2str(tmp{1,j})]) = tmp{i,j};
                end
            end
        end
    end
    clearvars tmp
end
if isempty(LAGS)
    % set default values from EGI
    LAGS.NA300.S1000 = 8;
    LAGS.NA300.S500 = 18;
    LAGS.NA300.S250 = 36;
    LAGS.NA300.S125 = 72;
    LAGS.NA400.S1000 = 36;
    LAGS.NA400.S500 = 66;
    LAGS.NA400.S250 = 112;
    LAGS.NA405.S1000 = 36;
    LAGS.NA405.S500 = 66;
    LAGS.NA405.S250 = 112;
    LAGS.NA410.S1000 = 13;
    LAGS.NA410.S500 = 34;
    LAGS.NA410.S250 = 76;
    LAGS.NONE.S1000 = 0;
    LAGS.NONE.S500 = 0;
    LAGS.NONE.S250 = 0;
    LAGS.NONE.S125 = 0;
end

if isfield(LAGS,amplifier) & isfield(LAGS.(amplifier),['S' num2str(samplingrate)])
    lag = LAGS.(amplifier).(['S' num2str(samplingrate)]);
    cfg.EXTRA.LAGS = LAGS;
    return
end
if isfield(cfg,'MAIN') & isfield(cfg.MAIN,'auto') & cfg.MAIN.auto == 1
    lag = 0;
    return
end
if ~isfield(LAGS,amplifier)
    Prompt = {'Amplifier','amplifier'};
    Formats.type = 'list';
    Formats.style = 'popupmenu';
    Formats.format = 'input';
    Formats.items = cat(1,{'OTHER'},fieldnames(LAGS));
    [Answer,Cancelled] = inputsdlg(Prompt,'Select Amplifier',Formats);
    pause(0.2);
    if Cancelled == 0 & ~strcmp(Answer.amplifier,'OTHER')
        LAGS.(amplifier) = LAGS.(Answer.amplifier);
    else
        LAGS.(amplifier) = [];
    end
    clearvars Prompt Formats
end
if isfield(LAGS.(amplifier),['S' num2str(samplingrate)])
    lag = LAGS.(amplifier).(['S' num2str(samplingrate)]);
    cfg.EXTRA.LAGS = LAGS;
    if isfield(LAGS,'DefaultEmpty')
        LAGS = rmfield(LAGS,'DefaultEmpty');
    end
    save(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'),'LAGS');
    return
end

Prompt = {['Lag for ' num2str(samplingrate) ' samplingrate'], 'lag'};
Formats.type = 'edit';
Formats.format = 'integer';
Formats.limits = [0 999999999];
Formats.size = 80;
[Answer,Cancelled] = inputsdlg(Prompt,'Select Lag',Formats);
if Cancelled == 1 | isempty(Answer.lag)
    lag = 0;
    return
else
    pause(0.2);
    lag = Answer.lag;
    LAGS.(amplifier).(['S' num2str(samplingrate)]) = Answer.lag;
    cfg.EXTRA.LAGS = LAGS;
    if isfield(LAGS,'DefaultEmpty')
        LAGS = rmfield(LAGS,'DefaultEmpty');
    end
    save(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'),'LAGS');
end

