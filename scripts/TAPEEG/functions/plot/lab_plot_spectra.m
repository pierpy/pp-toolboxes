% Helper function to plot spectras
%
% lab_plot_spectra(PLOT)
%
% written by F. Hatz 2012

function lab_plot_spectra(PLOT)

if ~exist('PLOT','var')
    [PLOT,skipprocessing] = lab_set_plot_spectra;
    if skipprocessing == 1
        return
    else
        pause(0.2);
    end
    figure1 = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none');
    if ~isempty(PLOT.epochsnumber)
        set(figure1,'Name',[PLOT.patient '_E' num2str(PLOT.epochsnumber) '_Channel ' num2str(PLOT.channels)]);
    else
        set(figure1,'Name',[PLOT.patient '_' PLOT.epochs '_Channel ' num2str(PLOT.channels)]);
    end
    m1 = uimenu(figure1,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback','close;');
    m2 = uimenu(figure1,'Label','Edit');
    uimenu(m2,'Label','Average method','Callback',@(~,~)select_avg);
    uimenu(m2,'Label','Set limits','Callback',@(~,~)set_limits);
    uimenu(m2,'Label','Select epoch','Callback',@(~,~)select_epoch);
    uimenu(m2,'Label','Select channels','Callback',@(~,~)select_channels);
end

T = PLOT.T;
R = PLOT.R;
if isfield(PLOT,'patient')
    patient = PLOT.patient;
else
    patient = '';
end
if isfield(T,'Valid')
    Valid = T.Valid;
else
    Valid = '';
end

SpectAllM = T.Spect;

do_plot;

if isfield(PLOT,'qualityplot') & PLOT.qualityplot == true
    plot_quality(T,R,SpectAllM,Valid,patient);
elseif isfield(PLOT,'correctpf') & PLOT.correctpf == true
    plot_pf(T,SpectAllM,patient);
else
    title(['Background Activity ' regexprep(patient,'_',' ')]);
end

    function do_plot
        area(SpectAllM);
        set(gca,'XTick',1:length(PLOT.freqlabel));
        set(gca,'XTickLabel',PLOT.freqlabel,'FontName','Times','fontsize',9);
        if isfield(PLOT,'Limits') & ~isempty(PLOT.Limits) & ...
                isfield(PLOT.Limits,'y1') & ~isempty(PLOT.Limits.y1) & ...
                isfield(PLOT.Limits,'y2') & ~isempty(PLOT.Limits.y2)
            set(gca,'YLim',[PLOT.Limits.y1 PLOT.Limits.y2]);
        end
    end
    function set_limits
        Prompt = {'Limit low','y1';'Limit high','y2'};
        Formats.type = 'edit';
        Formats.format = 'float';
        Formats.limits = [0 inf];
        if ~isfield(PLOT,'Limits') | ~isfield(PLOT.Limits,'y1') | isempty(PLOT.Limits.y1)
            PLOT.Limits.y1 = round(min(SpectAllM)*100)/100;
            PLOT.Limits.y2 = round(max(SpectAllM)*100)/100;
        end
        PLOT.Limits = inputsdlg(Prompt,'Channels',[Formats Formats],PLOT.Limits);
        clearvars Prompt Formats
        do_plot;
    end
    function select_epoch
       if strcmp(PLOT.epochs,'mean') | strcmp(PLOT.epochs,'median')
           PLOT.epochsnumber = 1:size(PLOT.Spect.SpectAll,1);
       end
       PLOT.epochs = 'single';
       Formats.type = 'list';
       Formats.style = 'listbox';
       Formats.limits = [0 2]; % multi-select
       Formats.size = [100 200];
       Formats.items = cellstr(num2str((1:size(PLOT.Spect.SpectAll,1))'));
       PLOT = inputsdlg({'Epoch number','epochsnumber'},'Epochs number',Formats,PLOT);
       clearvars Formats
       Spect = permute(PLOT.Spect.SpectAll(PLOT.epochsnumber,:,:),[3 2 1]);
       if strcmp(PLOT.source,'median')
           SpectAllM = median(median(Spect(PLOT.channels,:,:),1),3);
       else
           SpectAllM = mean(mean(Spect(PLOT.channels,:,:),1),3);
       end
       do_plot;
       set(figure1,'Name',[PLOT.patient '_E' num2str(PLOT.epochsnumber) '_Channel ' num2str(PLOT.channels)]);
   end
   function select_channels
       Formats.type = 'list';
       Formats.style = 'listbox';
       Formats.limits = [0 2]; % multi-select
       Formats.size = [100 200];
       if isfield(PLOT.Spect,'channels')
           Formats.items = PLOT.Spect.channels;
       elseif isfield(PLOT.Spect,'names')
           Formats.items = PLOT.Spect.names;
       else
           Formats.items = cellstr(num2str((1:size(PLOT.Spect.SpectAll,3))'));
       end
       PLOT = inputsdlg({'Channels','channels'},'Channels',Formats,PLOT);
       clearvars Formats
       if strcmp(PLOT.epochs,'mean')
           Spect = PLOT.Spect.SpectAllMean;
       elseif strcmp(PLOT.epochs,'median')
           Spect = PLOT.Spect.SpectAllMedian;
       else
           Spect = permute(PLOT.Spect.SpectAll(PLOT.epochsnumber,:,:),[3 2 1]);
       end
       if strcmp(PLOT.source,'median')
           SpectAllM = median(median(Spect(PLOT.channels,:,:),1),3);
       else
           SpectAllM = mean(mean(Spect(PLOT.channels,:,:),1),3);
       end
       do_plot;
       if ~isempty(PLOT.epochsnumber)
           set(figure1,'Name',[PLOT.patient '_E' num2str(PLOT.epochsnumber) '_Channel ' num2str(PLOT.channels)]);
       else
           set(figure1,'Name',[PLOT.patient '_' PLOT.epochs '_Channel ' num2str(PLOT.channels)]);
       end
   end
   function select_avg
       Formats.type = 'list';
       Formats.style = 'popupmenu';
       Formats.format = 'input';
       Formats.items = {'mean','median'};
       PLOT = inputsdlg({'Average method','source'},'Channels',Formats,PLOT);
       if strcmp(PLOT.epochs,'mean')
           Spect = PLOT.Spect.SpectAllMean;
       elseif strcmp(PLOT.epochs,'median')
           Spect = PLOT.Spect.SpectAllMedian;
       else
           Spect = permute(PLOT.Spect.SpectAll(PLOT.epochsnumber,:,:),[3 2 1]);
       end
       if strcmp(PLOT.source,'median')
           SpectAllM = median(median(Spect(PLOT.channels,:,:),1),3);
       else
           SpectAllM = mean(mean(Spect(PLOT.channels,:,:),1),3);
       end
   end

end

function plot_quality(T,R,SpectAllM,Valid,patient)
   SpectAllF = T.SpectF;
   if ~isempty(Valid)
       title(['Background Activity ' regexprep(patient,'_',' ') ' (' Valid ')']);
   else
       title(['Background Activity ' regexprep(patient,'_',' ')]);
   end
   hold on
   Lsize = max(SpectAllM);
   freqlow = T.areapower(1,1);
   freqhigh = T.areapower(1,2);
   peakfreq = T.peakfreq(1,1);
   if peakfreq > 0
       Values{1,1} = plot([peakfreq;peakfreq],[0;Lsize*1.1],'r');
       Values{2,1} = num2str(R.peakfreq(1),'%04.2f');
   else
       Values{1,1} = plot([freqlow;freqlow],[0;Lsize*0.05],'r');
       Values{2,1} = 'NaN';
   end
   Values{1,2} = plot([freqlow;freqlow],[0;Lsize*0.7],'--r');
   Values{2,2} = num2str(R.areapower(1),'%04.2f');
   plot([freqhigh;freqhigh],[0;Lsize*0.7],'--r');
   
   freqlow = T.cogareapower(1,1);
   freqhigh = T.cogareapower(1,2);
   cogfreq = T.cogfreq(1,1);
   if ~isnan(cogfreq)
       Values{1,3} = plot([cogfreq;cogfreq],[0;Lsize],'g');
       Values{2,3} = num2str(R.cogfreq(1),'%04.2f');
   else
       Values{1,3} = plot([freqlow;freqlow],[0;Lsize*0.2],'g');
       Values{2,3} = 'NaN';
   end
   Values{1,4} = plot([freqlow;freqlow],[0;Lsize*0.5],'--g');
   Values{2,4} = num2str(R.cogareapower(1),'%04.2f');
   plot([freqhigh;freqhigh],[0;Lsize*0.5],'--g');
   
   if ~strcmp(Values{2,1},'NaN')
       minfreqI = T.peak2min(1,1);
       minfreq = T.peak2min(1,2);
       clearvars tmp
       Values{1,5} = plot([minfreqI(1);minfreqI(1)],[0;minfreq(1)],'c');
       Values{2,5} = num2str(R.peak2min(1),'%04.2f');
   else
       Values{1,5} = plot([1;1],[0;Lsize*0.2],'g');
       Values{2,5} = 'NaN';
   end
   if isfield(R,'bplast') & R.bplast(1) > 0
       freqlow = T.bplast(1,1);
       freqhigh = T.bplast(1,2);
       Values{1,6} = plot([freqlow;freqlow],[0;Lsize*0.1],'y');
       Values{2,6} = num2str(R.bplast(1),'%04.2f');
       plot([freqhigh;freqhigh],[0;Lsize*0.1],'y');
   end
   if size(Values,2) == 5
       legend([Values{1,1} Values{1,2} Values{1,3} Values{1,4} Values{1,5}], ...
           ['Peak Frequency: ' Values{2,1} ' Hz ('  num2str(R.peakamp(1),'%04.2f') ')'], ...
           ['BandPower (PF): ' Values{2,2}], ...
           ['Median Frequency: ' Values{2,3} ' Hz'], ...
           ['BandPower (MF): ' Values{2,4}], ...
           ['Ratio Peak2Min: ' Values{2,5}]);
   elseif size(Values,2) == 6
       legend([Values{1,1} Values{1,2} Values{1,3} Values{1,4} Values{1,5} Values{1,6}], ...
           ['Peak Frequency: ' Values{2,1} ' Hz ('  num2str(R.peakamp(1),'%04.2f') ')'], ...
           ['BandPower (PF): ' Values{2,2}], ...
           ['Median Frequency: ' Values{2,3} ' Hz'], ...
           ['BandPower (MF): ' Values{2,4}], ...
           ['Ratio Peak2Min: ' Values{2,5}], ...
           ['BandPower ' num2str(SpectAllF(1,freqlow)) '-' ...
           num2str(SpectAllF(1,freqhigh)) 'Hz: ' Values{2,6}]);
   end
   clearvars legends Values freqlow freqhigh
end

function plot_pf(T,SpectAllM,patient)
   title(['Background Activity ' regexprep(patient,'_',' ')]);
   hold on
   Lsize = max(SpectAllM);
   freqlow = T.areapower(1,1);
   freqhigh = T.areapower(1,2);
   peakfreq = T.peakfreq(1,1);
   if peakfreq > 0
       plot([peakfreq;peakfreq],[0;Lsize*1.1],'r');
   end
   plot([freqlow;freqlow],[0;Lsize*0.7],'--r');
   plot([freqhigh;freqhigh],[0;Lsize*0.7],'--r');
   clearvars freqlow freqhigh
end