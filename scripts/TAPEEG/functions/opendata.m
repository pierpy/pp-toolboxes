% Script to open data in different formats
% 
% The script asks for a file to open and loads data by function lab_read_* 
% whereas * is replaced by file extension. Supported file formats are
% defined by files located in folder 'iofile'
%
% Written by F.Hatz, Neurology Basel 2012

[data,header,cfg] = lab_read_data;
if isempty(cfg)
    clearvars cfg
end
if isempty(header)
    clearvars header
end
if isempty(data) 
    clearvars data
    if exist('cfg','var') & ischar(cfg) & length(cfg) > 4 & strcmp(cfg(end-3:end),'.mat')
        load(cfg);
    end
end
if exist('cfg','var') & ischar(cfg)
    [~,~,Format] = lab_filename(cfg);
    if strcmp(Format,'xlsx') | strcmp(Format,'xls')
        selectmode = questdlg('Read Excel-File as statistics','Read Excel','Yes','No','No');
        if strcmp(selectmode,'Yes')
            [data,header,result,factors,cfg] = lab_read_statistics(data,1,1,1,0,1);
        end
        clearvars selectmode
    else
        clearvars cfg
    end
    clearvars Format
end
if exist('data','var') & ~isempty(data) & size(data,3) == 3 & isa(data,'uint8')
    figure('Color',[0 0 0],'MenuBar','None','Name','Image','NumberTitle','off');
    ScreenSize = get(0,'ScreenSize');
    ScreenSize(1) = ScreenSize(3) * 0.05;
    ScreenSize(2) = ScreenSize(4) * 0.95;
    ScreenSize(3) = ScreenSize(3) * 0.9;
    ScreenSize(4) = ScreenSize(4) * 0.9;
    if size(data,1) < ScreenSize(4) & size(data,2) < ScreenSize(3)
        ScreenSize(3) = size(data,2);
        ScreenSize(4) = size(data,1);
    end
    ScreenSize(2) = ScreenSize(2) - ScreenSize(4);
    set(gcf,'position',ScreenSize);
    clearvars ScreenSize
    image(data);
    set(gca, 'position', [0 0 1 1], 'visible', 'off');
    set(gcf,'PaperPositionMode','auto');
    axis off
    axis image
end



