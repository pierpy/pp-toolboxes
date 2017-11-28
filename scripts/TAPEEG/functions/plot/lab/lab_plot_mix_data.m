function [DATAout,PLOTout] = lab_plot_mix_data(DATA,PLOT,PlotList)
    
DATAout = [];
PLOTout = [];
for plotnr = 1:size(PlotList,2)
    plotnr2 = PlotList{plotnr};
    for i = plotnr2
        Mode = PLOT(i).Mode;
        if length(DATAout) < plotnr | ~isfield(DATAout(plotnr),Mode) | ...
                ~isfield(DATAout(plotnr).(Mode),'Nr') | isempty(DATAout(plotnr).(Mode).Nr)
            DATAout(plotnr,1).(Mode).Nr = 0; %#ok<AGROW>
        end
        if DATAout(plotnr,1).(Mode).Nr == 0
            DATAout(plotnr,1).(Mode).data = DATA(i).data; %#ok<AGROW>
            DATAout(plotnr,1).(Mode).Nr = 1; %#ok<AGROW>
            DATAout(plotnr,1).(Mode).measure = DATA(i).measure; %#ok<AGROW>
            DATAout(plotnr,1).(Mode).subject = DATA(i).subject; %#ok<AGROW>
            PLOTout(plotnr,1).(Mode).MinValue = PLOT(i).MinValue; %#ok<AGROW>
            PLOTout(plotnr,1).(Mode).MaxValue = PLOT(i).MaxValue; %#ok<AGROW>
            if isfield(PLOT,'Legend') & ~isempty(PLOT(i).Legend)
                PLOTout(plotnr,1).(Mode).Legend = PLOT(i).Legend; %#ok<AGROW>
            else
                PLOTout(plotnr,1).(Mode).Legend.Color = PLOT(i).Color; %#ok<AGROW>
                PLOTout(plotnr,1).(Mode).Legend.Text = {num2str(PLOT(i).MinValue,3);num2str(PLOT(i).MaxValue,3)}; %#ok<AGROW>
                PLOTout(plotnr,1).(Mode).Legend.Mode = 'Colorbar'; %#ok<AGROW>
            end
            datatmp = DATA(i).data;
        else
            if DATAout(plotnr,1).(Mode).Nr == 1
                DATAout(plotnr,1).(Mode).data = (DATAout(plotnr,1).(Mode).data - PLOTout(plotnr,1).(Mode).MinValue) / ...
                    (PLOTout(plotnr,1).(Mode).MaxValue - PLOTout(plotnr,1).(Mode).MinValue); %#ok<AGROW>
                DATAout(plotnr,1).(Mode).data(DATAout(plotnr,1).(Mode).data>1) = 1; %#ok<AGROW>
                DATAout(plotnr,1).(Mode).data(DATAout(plotnr,1).(Mode).data<0) = 0; %#ok<AGROW>
            end
            datatmp = (DATA(i).data - PLOT(i).MinValue) / (PLOT(i).MaxValue - PLOT(i).MinValue);
            datatmp(datatmp>1) = 1;
            datatmp(datatmp<0) = 0;
            datatmp = datatmp + DATAout(plotnr,1).(Mode).Nr;
            DATAout(plotnr,1).(Mode).data(datatmp>DATAout(plotnr,1).(Mode).Nr) = datatmp(datatmp>DATAout(plotnr,1).(Mode).Nr); %#ok<AGROW>
            if ~strcmp(DATAout(plotnr,1).(Mode).measure,DATA(i).measure)
                DATAout(plotnr,1).(Mode).measure = [DATAout(plotnr,1).(Mode).measure '_' DATA(i).measure]; %#ok<AGROW>
            end
            if ~strcmp(DATAout(plotnr,1).(Mode).subject,DATA(i).subject)
                DATAout(plotnr,1).(Mode).subject = [DATAout(plotnr,1).(Mode).subject '_' DATA(i).subject]; %#ok<AGROW>
            end
            DATAout(plotnr,1).(Mode).Nr = DATAout(plotnr,1).(Mode).Nr + 1; %#ok<AGROW>
            PLOTout(plotnr,1).(Mode).MinValue = 0; %#ok<AGROW>
            PLOTout(plotnr,1).(Mode).MaxValue = DATAout(plotnr,1).(Mode).Nr; %#ok<AGROW>
            if isfield(PLOT,'Legend') & ~isempty(PLOT(i).Legend)
                Legend = PLOT(i).Legend;
            else
                Legend.Color = PLOT(i).Color;
                Legend.Text = {num2str(PLOT(i).MinValue,3);num2str(PLOT(i).MaxValue,3)};
                Legend.Mode = 'Colorbar';
            end
            PLOTout(plotnr,1).(Mode).Legend = cat(1,PLOTout(plotnr,1).(Mode).Legend,Legend); %#ok<AGROW>
        end
        if isfield(PLOT,'Size') & ~isempty(PLOT(i).Size)
            if ~isfield(PLOTout(plotnr,1).(Mode),'Size')
                PLOTout(plotnr,1).(Mode).Size = ones(1,length(DATAout(plotnr,1).(Mode).data)); %#ok<AGROW>
                PLOTout(plotnr,1).(Mode).SizeE = ones(1,length(DATAout(plotnr,1).(Mode).data)); %#ok<AGROW>
            end
            if strcmp(Mode,'Connections')
                PLOTout(plotnr,1).(Mode).Size(DATA(i).data>0) = 1; %#ok<AGROW>
                PLOTout(plotnr,1).(Mode).SizeE(DATA(i).data>0) = PLOT(i).Size; %#ok<AGROW>
            else
                PLOTout(plotnr,1).(Mode).Size(DATA(i).data>0) = PLOT(i).Size; %#ok<AGROW>
            end
        end
        if isfield(PLOT,'Alpha') & ~isempty(PLOT(i).Alpha)
            if ~isfield(PLOTout(plotnr,1).(Mode),'Alpha')
                PLOTout(plotnr,1).(Mode).Alpha = ones(1,length(DATAout(plotnr,1).(Mode).data)); %#ok<AGROW>
            end
            PLOTout(plotnr,1).(Mode).Alpha(DATA(i).data>0) = PLOT(i).Alpha; %#ok<AGROW>
        end
        if ~isfield(PLOTout(plotnr,1).(Mode),'Color') | isempty(PLOTout(plotnr,1).(Mode).Color)
            PLOTout(plotnr,1).(Mode).Color = PLOT(i).Color; %#ok<AGROW>
        else
            PLOTout(plotnr,1).(Mode).Color = cat(1,PLOTout(plotnr,1).(Mode).Color,PLOT(i).Color); %#ok<AGROW>
        end
    end
end

end