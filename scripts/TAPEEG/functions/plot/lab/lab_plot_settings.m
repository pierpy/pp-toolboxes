function PLOT = lab_plot_settings(PLOT,DATA,doIS)

if ~exist('doIS','var')
    doIS = 0;
end
if ~exist('DATA','var')
    DATA = [];
end
if ~exist('PLOT','var')
    PLOT = [];
end
if isempty(PLOT) & isempty(DATA)
    return
end
doconnections = false;
if isfield(DATA,'data') & ~isempty(DATA(1).data)
    YData = size(DATA,1);
    maxsubject = size(DATA,2);
    StatFlag = zeros(1,maxsubject);
    if isfield(DATA,'subject') & ~isempty(DATA(1).subject)
        for datanr = 1:maxsubject
            if strcmp(DATA(1,datanr).subject,'p') | strcmp(DATA(1,datanr).subject,'mxp') | ...
                    strcmp(DATA(1,datanr).subject,'mxpV') | strcmp(DATA(1,datanr).subject,'mxpV2') | ...
                    strcmp(DATA(1,datanr).subject,'mxpC') | strcmp(DATA(1,datanr).subject,'mxpS')
                StatFlag(datanr) = 1;
                %             elseif strfind(DATA(1,datanr).subject,'_p')
                %                 StatFlag(datanr) = 1;
                %             elseif strfind(DATA(1,datanr).subject,'_mxp')
                %                 StatFlag(datanr) = 1;
                %             elseif strfind(DATA(1,datanr).subject,'_mxpV')
                %                 StatFlag(datanr) = 1;
                %             elseif strfind(DATA(1,datanr).subject,'_mxpV2')
                %                 StatFlag(datanr) = 1;
                %             elseif strfind(DATA(1,datanr).subject,'_mxpC')
                %                 StatFlag(datanr) = 1;
                %             elseif strfind(DATA(1,datanr).subject,'_mxpS')
                %                 StatFlag(datanr) = 1;
            elseif (strcmp(DATA(1,datanr).subject(end),'T') | strcmp(DATA(1,datanr).subject(end),'F') |  ...
                    strcmp(DATA(1,datanr).subject(end),'Z'))
                StatFlag(datanr) = 2;
            elseif strcmp(DATA(1,datanr).subject(end),'R')
                StatFlag(datanr) = 3;
            end
        end
    end
    tmp = find(StatFlag == 1,1,'first');
    if ~isempty(tmp)
        PosFlag = zeros(size(DATA));
        PosFlag(:,1:tmp) = 1;
        StatFlag = repmat(StatFlag,size(DATA,1),1);
        DATA = DATA';
        StatFlag = StatFlag';
        PosFlag = PosFlag';
        DATA = DATA(:);
        StatFlag = StatFlag(:);
        PosFlag = PosFlag(:);
        Mode = zeros(size(DATA,1),1);
        for datanr = 1:size(DATA,1)
            if isfield(DATA(datanr,1),'nodes') & DATA(datanr,1).nodes == true
                Mode(datanr) = 1;
            end
            if isfield(DATA(datanr,1),'connections') & DATA(datanr,1).connections == true
                Mode(datanr) = 2;
            end
        end
        for datanr = 2:size(DATA,1)
            if Mode(datanr) == 0 & ~strcmp(DATA(datanr).subject,DATA(datanr-1).subject)
                PosFlag(datanr,1) = 0;
            elseif Mode(datanr) == 1 & Mode(datanr-1) == 2
                PosFlag(datanr,1) = 0;
            end
        end
    else
        PosFlag = ones(size(DATA,1),1);
        Mode = zeros(size(DATA,1),1);
        for datanr = 1:size(DATA,1)
            if isfield(DATA(datanr,1),'nodes') & DATA(datanr,1).nodes == true
                Mode(datanr) = 1;
            end
            if isfield(DATA(datanr,1),'connections') & DATA(datanr,1).connections == true
                Mode(datanr) = 2;
            end
        end
        for datanr = 2:size(DATA,1)
            if Mode(datanr) == 0 & ~strcmp(DATA(datanr).subject,DATA(datanr-1).subject)
                PosFlag(datanr,1) = 0;
            elseif Mode(datanr) == 1 & Mode(datanr-1) == 2
                PosFlag(datanr,1) = 0;
            end
        end
        StatFlag = zeros(size(DATA,1),1);
        DATA = DATA(:,1);
    end
    clearvars tmp
end
if isempty(PLOT) | (~isempty(DATA) & length(PLOT) ~= length(DATA))
    maxdata = length(DATA);
    PLOT = [];
    for datanr = 1:maxdata
        % set default settings for every datarow
        if ~isfield(DATA(1),'labels') | isempty(DATA(1).labels) | DATA(1).labels == false
            if StatFlag(datanr) > 0
                PLOT(datanr,1).Name = regexprep(DATA(datanr).name,'_',' '); %#ok<AGROW>
            else
                PLOT(datanr,1).Name = regexprep(DATA(datanr).measure,'_',' '); %#ok<AGROW>
            end
            if StatFlag(datanr) == 1
                PLOT(datanr,1).Inverse = true; %#ok<AGROW>
                PLOT(datanr,1).MinValue = 0; %#ok<AGROW>
                PLOT(datanr,1).MaxValue = 0; %#ok<AGROW>
                PLOT(datanr,1).Threshold1 = 0.01; %#ok<AGROW>
                PLOT(datanr,1).Threshold2 = 0.05; %#ok<AGROW>
                PLOT(datanr,1).ColorMode = 'color'; %#ok<AGROW>
                PLOT(datanr,1).Color1 = [1 0 0]; %#ok<AGROW>
                PLOT(datanr,1).Color2 = [1 0 1]; %#ok<AGROW>
                if PosFlag(datanr) == 0 & datanr > 1 & Mode(datanr) == Mode(datanr-1)
                    PLOT(datanr-1,1).Color1 = [0 1 0]; %#ok<AGROW>
                    PLOT(datanr-1,1).Color2 = [0 1 1]; %#ok<AGROW>
                end
                if PosFlag(datanr) == 0
                    PLOT(datanr,1).AddPlot = true; %#ok<AGROW>
                else
                    PLOT(datanr,1).AddPlot = false; %#ok<AGROW>
                end
                if Mode(datanr) == 2
                    PLOT(datanr,1).Mode = 'Connections'; %#ok<AGROW>
                    PLOT(datanr,1).Size = 1; %#ok<AGROW>
                    if doIS > 0
                        PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
                    end
                elseif doIS == false | Mode(datanr) == 1
                    PLOT(datanr,1).Mode = 'Nodes'; %#ok<AGROW>
                    PLOT(datanr,1).Size = 1; %#ok<AGROW>
                    if doIS > 0
                        PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
                    end
                else
                    PLOT(datanr,1).Mode = 'Surface'; %#ok<AGROW>
                    if doIS > 0
                        PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
                        PLOT(datanr,1).Size = []; %#ok<AGROW>
                    else
                        PLOT(datanr,1).Size = 1; %#ok<AGROW>
                    end
                end
            else
                PLOT(datanr,1).Inverse = false; %#ok<AGROW>
                if StatFlag(datanr) == 2
                    PLOT(datanr,1).MinValue = -5; %#ok<AGROW>
                    PLOT(datanr,1).MaxValue = 5; %#ok<AGROW>
                    PLOT(datanr,1).Threshold1 = []; %#ok<AGROW>
                    PLOT(datanr,1).Threshold2 = []; %#ok<AGROW>
                    PLOT(datanr,1).ColorMode = 'bluered'; %#ok<AGROW>
                    PLOT(datanr,1).Color1 = [1 1 1]; %#ok<AGROW>
                elseif StatFlag(datanr) == 3
                    PLOT(datanr,1).MinValue = -1; %#ok<AGROW>
                    PLOT(datanr,1).MaxValue = 1; %#ok<AGROW>
                    PLOT(datanr,1).Threshold1 = []; %#ok<AGROW>
                    PLOT(datanr,1).Threshold2 = []; %#ok<AGROW>
                    PLOT(datanr,1).ColorMode = 'bluered'; %#ok<AGROW>
                    PLOT(datanr,1).Color1 = [1 1 1]; %#ok<AGROW>
                else
                    if min(DATA(datanr).data) < 0
                        if min(DATA(datanr).data) > -1
                            PLOT(datanr,1).MinValue = -max(max(abs(DATA(datanr).data))); %#ok<AGROW>
                        else
                            PLOT(datanr,1).MinValue = -ceil(max(max(abs(DATA(datanr).data)))); %#ok<AGROW>
                        end
                    else
                        PLOT(datanr,1).MinValue = min(min(DATA(datanr).data)); %#ok<AGROW>
                    end
                    if max(abs(DATA(datanr).data)) < 1
                        PLOT(datanr,1).MaxValue = max(max(abs(DATA(datanr).data))); %#ok<AGROW>
                    else
                        PLOT(datanr,1).MaxValue = ceil(max(max(abs(DATA(datanr).data)))); %#ok<AGROW>
                    end
                    PLOT(datanr,1).Threshold1 = []; %#ok<AGROW>
                    PLOT(datanr,1).Threshold2 = []; %#ok<AGROW>
                    if min(DATA(datanr).data) < 0
                        PLOT(datanr,1).ColorMode = 'bluered'; %#ok<AGROW>
                        PLOT(datanr,1).Color1 = [1 1 1]; %#ok<AGROW>
                    else
                        PLOT(datanr,1).ColorMode = 'color'; %#ok<AGROW>
                        PLOT(datanr,1).Color1 = [1 0 0]; %#ok<AGROW>
                    end
                end
                PLOT(datanr,1).Color2 = [1 1 1]; %#ok<AGROW>
                if PosFlag(datanr) == 0
                    PLOT(datanr,1).AddPlot = true; %#ok<AGROW>
                else
                    PLOT(datanr,1).AddPlot = false; %#ok<AGROW>
                end
                if Mode(datanr) == 2
                    PLOT(datanr,1).Mode = 'Connections'; %#ok<AGROW>
                    if doIS > 0
                        PLOT(datanr,1).Color1 = [1 0 0]; %#ok<AGROW>
                        PLOT(datanr,1).Size = 1; %#ok<AGROW>
                        PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
                    else
                        PLOT(datanr,1).Color1 = [0 0 1]; %#ok<AGROW> %[0.1686    0.5059    0.3373];
                        PLOT(datanr,1).Size = 1; %#ok<AGROW>
                    end
                elseif Mode(datanr) == 1
                    PLOT(datanr,1).Mode = 'Nodes'; %#ok<AGROW>
                    if doIS > 0
                        PLOT(datanr,1).Size = 1; %#ok<AGROW>
                        PLOT(datanr,1).Color1 = [0 0 0]; %#ok<AGROW>
                        PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
                    else
                        PLOT(datanr,1).ColorMode = 'color'; %#ok<AGROW>
                        PLOT(datanr,1).Color1 = [1 0 0]; %#ok<AGROW>
                    end
                else
                    PLOT(datanr,1).Mode = 'Surface'; %#ok<AGROW>
                    if doIS == 0
                        PLOT(datanr,1).Size = 1; %#ok<AGROW>
                    else
                        PLOT(datanr,1).Size = []; %#ok<AGROW>
                        PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
                    end
                end
            end
        else
            PLOT(datanr,1).Name = regexprep(DATA(datanr).name,'_',' '); %#ok<AGROW>
            PLOT(datanr,1).Channels = find(DATA(datanr).data == 1); %#ok<AGROW>
            if datanr == 1
                PLOT(datanr,1).Color1 = [0 1 0]; %#ok<AGROW>
            elseif datanr == 2
                PLOT(datanr,1).Color1 = [0 1 1]; %#ok<AGROW>
            elseif datanr == 3
                PLOT(datanr,1).Color1 = [1 0 1]; %#ok<AGROW>
            elseif datanr == 4
                PLOT(datanr,1).Color1 = [1 1 0]; %#ok<AGROW>
            else
                PLOT(datanr,1).Color1 = rand(1,3); %#ok<AGROW>
            end
            if doIS == 0
                PLOT(datanr,1).MaxDistance = 0; %#ok<AGROW>
                PLOT(datanr,1).WeightMaxDistance = 100; %#ok<AGROW>
            end
            if datanr > 1
                PLOT(datanr,1).AddPlot = true; %#ok<AGROW>
            else
                PLOT(datanr,1).AddPlot = false; %#ok<AGROW>
            end
            if doIS == 0
                PLOT(datanr,1).Mode = 'Nodes'; %#ok<AGROW>
                PLOT(datanr,1).Size = 1; %#ok<AGROW>
            else
                PLOT(datanr,1).Mode = 'Surface'; %#ok<AGROW>
                PLOT(datanr,1).Size = []; %#ok<AGROW>
                PLOT(datanr,1).Alpha = 1; %#ok<AGROW>
            end
        end
        if isfield(DATA(datanr),'connections') & DATA(datanr).connections == true
            PLOT(datanr,1).Mode = 'Connections'; %#ok<AGROW>
            doconnections = true;
        end
    end
    if doconnections == true
        fignr = [];
        for datanr = 1:length(PLOT)
            if PLOT(datanr,1).AddPlot == false
                fignr = datanr;
            end
            if strcmp(PLOT(datanr,1).Mode,'Surface') & (datanr - fignr) <= 1
                PLOT(datanr,1).Mode = 'Nodes'; %#ok<AGROW>
            end
        end
    end
else
    YData = size(PLOT,1);
    PLOT = PLOT';
    PLOT = PLOT(:);
end

FORMAT = {};
if isfield(PLOT,'Settings')
    tmp = [];
    tmp.type = 'check';
    tmp.callback = {@check_color,'Settings','Settings'};
    FORMAT{strcmp(fieldnames(PLOT),'Settings')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Name')
    FORMAT{strcmp(fieldnames(PLOT),'Name')} = 'text';
end
if isfield(PLOT,'Channels')
    FORMAT{strcmp(fieldnames(PLOT),'Channels')} = 'vector';
end
if isfield(PLOT,'ColorMode')
    tmp = [];
    tmp.type = 'list';
    tmp.style = 'popupmenu';
    tmp.format = 'input';
    tmp.items = {'color','bluered','autumn','bone','colorcube','cool','copper','gray', ...
        'hot','hsv','jet','iautumn','ibone','icolorcube','icool','icopper','igray','ihot','ihsv','ijet'};
    tmp.callback = {@check_colormode,'Color1','ColorMode','Color1'};
    FORMAT{strcmp(fieldnames(PLOT),'ColorMode')} = tmp;
    clearvars tmp
    if isfield(PLOT,'Color1')
        tmp = [];
        tmp.type = 'color';
        tmp.callback = {@check_color,'ColorMode','ColorMode','Color1'};
        FORMAT{strcmp(fieldnames(PLOT),'Color1')} = tmp;
        clearvars tmp
    end
end
if isfield(PLOT,'Mode')
    tmp = [];
    tmp.type = 'list';
    tmp.style = 'popupmenu';
    tmp.format = 'input';
    if doIS == 2
        tmp.items = {'Nodes','Surface','Volume'};
    else
        tmp.items = {'Nodes','Surface'};
    end
    if isfield(PLOT,'Size')
        tmp.callback = {@check_mode,'Size','Mode','Size',doIS};
    end
    FORMAT{strcmp(fieldnames(PLOT),'Mode')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Threshold1')
    tmp = [];
    tmp.type = 'edit';
    tmp.format = 'float';
    tmp.size = 50;
    FORMAT{strcmp(fieldnames(PLOT),'Threshold1')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Threshold2')
    tmp = [];
    tmp.type = 'edit';
    tmp.format = 'float';
    tmp.size = 50;
    tmp.callback = {@check_threshold2,'Threshold2','Threshold2','Threshold1','Inverse'};
    FORMAT{strcmp(fieldnames(PLOT),'Threshold2')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Color2')
    tmp = [];
    tmp.type = 'color';
    tmp.callback = {@check_color2,'Threshold2','Threshold2','Color2'};
    FORMAT{strcmp(fieldnames(PLOT),'Color2')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Size')
    tmp = [];
    tmp.type = 'edit';
    tmp.format = 'float';
    tmp.size = 50;
    FORMAT{strcmp(fieldnames(PLOT),'Size')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Alpha')
    tmp = [];
    tmp.type = 'edit';
    tmp.format = 'float';
    tmp.size = 50;
    tmp.limits = [0 1];
    FORMAT{strcmp(fieldnames(PLOT),'Alpha')} = tmp;
    clearvars tmp
end
if isfield(PLOT,'Valid')
    tmp = [];
    tmp.type = 'check';
    tmp.enable = 'inactive';
    FORMAT{strcmp(fieldnames(PLOT),'Valid')} = tmp;
    clearvars tmp
end

numdiags = ceil(size(PLOT,1) / 20);
for i = 1:numdiags
    Pstart = (i-1)*20+1;
    Pstop = i*20;
    if Pstop > size(PLOT,1)
        Pstop = size(PLOT,1);
    end
    if i == numdiags
        OPTIONS.ButtonNames = {'OK','Cancel'};
    else
        OPTIONS.ButtonNames = {'Next','Cancel'};
    end
    [answer,skipprocessing] = inputsdlg(PLOT(Pstart:Pstop),'Plot settings',FORMAT,[],OPTIONS);
    if skipprocessing == 1
        PLOT = [];
        return
    else
        PLOT(Pstart:Pstop) = answer;
        clearvars answer
    end
end
if size(PLOT,1) > YData
    PLOT = reshape(PLOT,size(PLOT,1)/YData,YData);
    PLOT = PLOT';
end

end

function Color1 = check_colormode(ColorMode,Color1)
    if ~strcmp(ColorMode,'color')
        Color1 = [1 1 1];
    end
end

function ColorMode = check_color(ColorMode,Color1)
    if min(Color1) ~= 1 | max(Color1) ~= 1
        ColorMode = 'color';
    end
end

function Threshold2 = check_threshold2(Threshold2,Threshold1,Inverse)
    if Inverse == true & ~isempty(Threshold2) & Threshold2 < Threshold1
        Threshold2 = [];
    elseif Inverse == false & ~isempty(Threshold2) & Threshold2 > Threshold1
        Threshold2 = [];
    end
end

function Threshold2 = check_color2(Threshold2,Color2)
    if min(Color2) == 1 & max(Color2) == 1
        Threshold2 = [];
    elseif isempty(Threshold2)
        Threshold2 = 1;
    end
end

function Size = check_mode(Mode,Size,doIS)
    if doIS > 0 & (strcmp(Mode,'Surface') | strcmp(Mode,'Volume'))
        Size = [];
    else
        Size = 1;
    end 
end