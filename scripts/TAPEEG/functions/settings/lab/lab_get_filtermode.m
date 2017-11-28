function settings = lab_get_filtermode(settings,Force)

if exist('Force','var') & Force == 1
    settings.filtorder = [];
end
if strcmp(settings.filtermode,'wavelet')
    settings.filtorder = [];
    if ~isfield(settings,'wavsmoothing') | isempty(settings.wavsmoothing)
        settings.wavsmoothing = 4;
    end
    Prompt(1,:) = {'Smoothing factor (wavelet)','wavsmoothing'};
    Formats(1,1).type = 'edit';
    Formats(1,1).format = 'float';
    Formats(1,1).limits = [0 inf];
    Formats(1,1).size = 30;
    [settings,Cancelled] = inputsdlg(Prompt,'Wavelet setting',Formats,settings);
    if Cancelled == 1
        settings.wavsmoothing = 4;
    end
elseif strcmp(settings.filtermode,'butter') | strcmp(settings.filtermode,'firls') | strcmp(settings.filtermode,'cheby')
    settings.wavsmoothing = [];
    if strcmp(settings.filtermode,'butter')
        if ~isfield(settings,'filtorder') | isempty(settings.filtorder) | settings.filtorder == 0
            settings.filtorder = 2;
        end
        Prompt(1,:) = {'Filter order','filtorder'};
    elseif strcmp(settings.filtermode,'firls')
        if ~isfield(settings,'filtorder') | isempty(settings.filtorder)
            settings.filtorder = 0;
        end
        Prompt(1,:) = {'Filter order (0=4.8*samplingrate)','filtorder'};
    elseif strcmp(settings.filtermode,'cheby')
        if ~isfield(settings,'filtorder') | isempty(settings.filtorder)
            settings.filtorder = 0;
        end
        Prompt(1,:) = {'Filter order (0=MinOrder)','filtorder'};
    end
    Formats(1,1).type = 'edit';
    Formats(1,1).format = 'float';
    Formats(1,1).limits = [0 inf];
    Formats(1,1).size = 30;
    [settings,Cancelled] = inputsdlg(Prompt,'Filter order',Formats,settings);
    if Cancelled == 1
        settings.filtorder = 0;
    end
else
    settings.filtorder = [];
    settings.wavsmoothing = [];
end