function lag = lab_define_egi_lag(samplingrate)
    
global Main_Path
if ~exist('samplingrate','var')
    samplingrate = [];
end
lag = [];

[~,Filepath] = uigetfile('*.bin','select signal1.bin in mff-folder');
if ~ischar(Filepath) | isempty(Filepath)
    return
end

InfoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info,'info.xml',Filepath);
amplifier = char(InfoObj.getAmpSerialNumber);
if isempty(amplifier)
    amplifier = 'DefaultEmpty';
end

if exist(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'),'file')
    load(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'))
else
    LAGS.NONE.S1000 = 0;
    LAGS.NONE.S500 = 0;
    LAGS.NONE.S250 = 0;
    LAGS.NONE.S125 = 0;
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
end

Prompt = {'Amplifier','amplifier'};
Formats.type = 'list';
Formats.style = 'popupmenu';
Formats.format = 'input';
Formats.items = fieldnames(LAGS);
[Answer,Cancelled] = inputsdlg(Prompt,'Select Amplifier',Formats);
pause(0.2);
if Cancelled == 1
    return
end
clearvars Prompt Formats

LAGS.(amplifier) = LAGS.(Answer.amplifier);
save(fullfile(fullfile(Main_Path,'.ignore'),'amplag.mat'),'LAGS');

if isfield(LAGS.(amplifier),['S' num2str(samplingrate)])
    lag = LAGS.(amplifier).(['S' num2str(samplingrate)]);
end

end