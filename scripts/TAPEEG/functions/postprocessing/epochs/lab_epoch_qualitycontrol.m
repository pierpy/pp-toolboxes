% Signal-space-plot of bandpower 8-12Hz (alpha), bad electrodes in every
% epoch are labeled green
%
% lab_epoch_qualitycontrol(epochs,header,filename)
%
% written by F. Hatz 2012

function lab_epoch_qualitycontrol(epochs,header,filename)

disp('    plot quality-check-plots');
[~,filepath,~,filenameS] = lab_filename(filename);

if ~isfield(header,'locs') | ~isfield(header.locs,'x')
    return
end
if isfield(header,'numdatachannels') & size(header.locs.x,2) ~= header.numdatachannels
    return
elseif ~isfield(header,'numdatachannels') & size(header.locs.x,2) > size(epochs.data,1)
    return
elseif ~isfield(header,'numdatachannels') & size(header.locs.x,2) <= size(epochs.data,1)
    header.numdatachannels = size(header.locs.x,2);
end

maxnumepochs = floor(size(epochs.markersvalid,2)/epochs.markersvalidshift);

settings = [];
for nepoch = 1:size(epochs.data,3)
    settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
    epochfft = [];
    for i = 1:size(epochs.data,1)
        [epochfft(i,:),~,freqs] = pmtm(detrend(epochs.data(i,:,nepoch)),4,size(epochs.data,2),header.samplingrate); %#ok<AGROW>
    end
    deltafreq = freqs(2,1)-freqs(1,1);
    lowfreq = round(8 / deltafreq)+1;
    highfreq = round(12 / deltafreq)+1;
    epochfft = epochfft(1:header.numdatachannels,lowfreq:highfreq);
    epochfft = max(epochfft,[],2);
    if ~isfield(settings,'LOCS')
        settings.LOCS = header.locs;
        settings.PLOT_file = [];
    end
    PLOT.MinValue = -max(abs(epochfft(:)));
    PLOT.MaxValue = max(abs(epochfft(:)));
    PLOT.Color = lab_create_cmap('r');
    settings = lab_plot_chans(epochfft,PLOT,settings);
    clearvars PLOT
    
    weights = epochs.badoutelectrodes(nepoch,:);
    weights(weights<=0) = -1.2;
    PLOT.MinValue = -1.2;
    PLOT.MaxValue = 1.2;
    PLOT.AddPlot = 1;
    PLOT.Color = lab_create_cmap('g');
    settings = lab_plot_chans(weights,PLOT,settings,1);
    clearvars PLOT
    
    text(-0.5,-0.5,['%good: ' num2str(epochs.goodsum(1,nepoch))],'HorizontalAlignment','left','VerticalAlignment','bottom');
    text(0,0.5,['Distribution 8 - 12 Hz (epoch ' num2str(nepoch) '/' num2str(size(epochs.data,3)) ' of ' num2str(maxnumepochs) ')'],'HorizontalAlignment','center','VerticalAlignment','top');
    set(gcf, 'Name',fullfile(filepath,[filenameS  '_E' num2str(nepoch)]),'NumberTitle','off');
    
    lab_print_figure(fullfile(filepath,[filenameS  '_E' num2str(nepoch) '.jpg']),settings.handleF);
    close(settings.handleF);
end