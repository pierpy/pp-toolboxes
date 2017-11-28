function [settings,skipprocessing] = lab_set_plot_spectra(settings)

disp ('   Ask for collect spectraldata settings')

skipprocessing = 0;

if ~exist('settings','var') | ~isfield(settings,'source')
    settings.file = '';
    settings.Spect = [];
    settings.Locs = [];
    settings.channels = [];
    settings.source = 'median';
    settings.epochs = 'median';
    settings.epochsnumber = [];
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Spectra-File','file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.mat;*.xls;*.xlsx','Spectra-File (*.mat,*.xls,*.xlsx)'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = 300;
Formats(end,1).callback = {@load_spectra,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'Spectra','Spect'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_spectra,'@ALL','@ALL'};

Prompt(end+1,:) = {'LOCS','Locs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_locs,'Locs','Locs','Spect'};

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Channels','channels'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [1 999999];
Formats(end,1).enable = 'inactive';
Formats(end,1).size = 300;
Formats(end,1).callback = {@set_channels,'@ALL','@ALL'};
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'','source'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'mean','median'};

Prompt(end+1,:) = {'',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 4];

Prompt(end+1,:) = {'Epochs','epochs'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'mean','median','single'};
Formats(end,1).callback = {@set_epochs,'@ALL','@ALL'};

Prompt(end+1,:) = {'Epoch number','epochsnumber'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [1 9999];
Formats(end,1).size = 30;

[settings,Cancelled] = inputsdlg(Prompt,'Plot spectrum',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
else
    if isfield(settings.Spect,'header') & isfield(settings.Spect.header,'patient')
        settings.patient = settings.Spect.header.patient;
    else
        settings.patient = lab_subjectname(settings.file);
    end
    if iscell(settings.Spect.SpectAllF)
        settings.freqlabel = settings.Spect.SpectAllF;
    else
        settings.freqlabel = cell(1,length(settings.Spect.SpectAllF));
        for i = 1:length(settings.Spect.SpectAllF)
            if round(settings.Spect.SpectAllF(i)) == settings.Spect.SpectAllF(i)
                settings.freqlabel(i) = num2cell(settings.Spect.SpectAllF(i));
            end
        end
        tmp = find(~cellfun(@isempty,settings.freqlabel));
        if ~isempty(tmp) & length(tmp) >25
            for i = 2:2:length(tmp)
                settings.freqlabel{1,tmp(i)} = [];
            end
        elseif isempty(tmp) | length(tmp) < 5
            settings.freqlabel = cell(1,length(settings.Spect.SpectAllF));
            for i = 1:floor(length(settings.Spect.SpectAllF)/25):length(settings.Spect.SpectAllF)
                settings.freqlabel{1,i} = settings.Spect.SpectAllF(i);
            end
        end
        clearvars tmp
    end
    settings.R = [];
    if strcmp(settings.epochs,'mean')
        Spect = settings.Spect.SpectAllMean;
        settings.epochsnumber = 1:size(settings.Spect.SpectAll,1);
    elseif strcmp(settings.epochs,'median')
        Spect = settings.Spect.SpectAllMedian;
        settings.epochsnumber = 1:size(settings.Spect.SpectAll,1);
    else
        if settings.epochsnumber > size(settings.Spect.SpectAll,1)
            settings.epochsnumber = size(settings.Spect.SpectAll,1);
        end
        if settings.epochsnumber < 1
            settings.epochsnumber = 1;
        end
        Spect = permute(settings.Spect.SpectAll(settings.epochsnumber,:,:),[3 2 1]);
    end
    if strcmp(settings.source,'median')
        settings.T.Spect = median(Spect(settings.channels,:),1);
    else
        settings.T.Spect = mean(Spect(settings.channels,:),1);
    end
end

end

function settings = load_spectra(settings)
   settings.Spect = [];
   if ~isempty(settings.file) & exist(settings.file,'file')
       [~,~,Format] = lab_filename(settings.file);
       switch Format
           case 'mat'
               try %#ok<TRYNC>
                   settings.Spect = load(settings.file);
                   if isfield(settngs.Spect,'header') & isfield(settings.Spect.header,'locs') & ~isempty(settings.Spect.header.locs)
                       settings.Spect.channels = settings.Spect.header.locs.labels';
                   else
                       for j = 1:size(settings.Spect.SpectAllMean,1)
                           settings.Spect.channels{j,1} = ['Channel' num2str(j)];
                       end
                   end
               end
           case {'xls','xlsx'}
               try %#ok<TRYNC>
                   result = lab_read_xls(settings.file);
                   ischannels = false;
                   if ischar(result{1,1}) & (~isnumeric(result{1,2}) | ...
                           strcmp(result{1,1},'Patient') | strcmp(result{1,1},'Region') | ...
                           strcmp(result{1,1},'Channel'))
                       if strcmp(result{1,1},'Channel')
                           ischannels = true;
                       end
                       if strcmp(result{1,1},'Patient')
                           settings.source = 'single';
                       end
                       if min(cellfun(@isnumeric,result(1,2:end))) == 1
                           freqlabel = cell2mat(result(1,2:end));
                       else
                           freqlabel = result(1,2:end);
                           for j = 1:size(freqlabel,2)
                               if isnumeric(freqlabel{1,j})
                                   freqlabel{1,j} = num2str(freqlabel{1,j});
                               end
                           end
                       end
                       channels = result(2:end,1);
                       for j = 1:size(channels,1)
                           if isnumeric(channels{j,1})
                               channels{j,1} = ['Spectrum' num2str(channels{j,1})];
                           end
                       end
                       result = cell2mat(result(2:end,2:end));
                   elseif ischar(result{1,1}) & isnumeric(result{1,2})
                       channels = result(1:end,1);
                       result = cell2mat(result(1:end,2:end));
                       settings2.lowfreq = 1;
                       settings2.freqbin = 0.25;
                       Formats.format = 'float';
                       Formats.size = 40;
                       Formats(2,1).format = 'float';
                       Formats(2,1).size = 40;
                       settings2 = inputsdlg({'Lowest Frequency','lowfreq';'Freqbin','freqbin'}, ...
                           'Frequency range',Formats,settings2,2);
                       freqlabel = settings2.lowfreq:settings2.freqbin:ceil(settings2.freqbin*size(result,2));
                       freqlabel = freqlabel(1,1:size(result,2));
                       clearvars settings2
                   else
                       result = cell2mat(result);
                       for j = 1:size(result,1)
                           channels{j,1} = ['Spectrum_' num2str(j)];
                       end
                       settings2.lowfreq = 1;
                       settings2.freqbin = 0.25;
                       Formats.format = 'float';
                       Formats.size = 40;
                       Formats(2,1).format = 'float';
                       Formats(2,1).size = 40;
                       settings2 = inputsdlg({'Lowest Frequency','lowfreq';'Freqbin','freqbin'}, ...
                           'Frequency range',Formats,settings2,2);
                       freqlabel = settings2.lowfreq:settings2.freqbin:ceil(settings2.freqbin*size(result,2));
                       freqlabel = freqlabel(1,1:size(result,2));
                       clearvars settings2
                   end
                   settings.Spect.SpectAllF = freqlabel;
                   settings.Spect.SpectAllMean = result;
                   settings.Spect.SpectAllMedian = result;
                   settings.Spect.SpectAll = permute(result,[3 2 1]);
                   if ischannels == true
                       settings.Spect.channels = channels;
                   else
                       settings.Spect.names = channels;
                   end
                   settings.channels = 1:length(channels);
               end
       end
   end
   if isfield(settings.Spect,'header') & isfield(settings.Spect.header,'locs') & ...
           ~isempty(settings.Spect.header.locs)
       settings.Locs = settings.Spect.header.locs;
   end
   if strcmp(settings.source,'single')
       settings.mappingsBA = 1;
   else
       if isfield(settings.Spect,'SpectAll')
           BA = lab_get_bactivity(size(settings.Spect.SpectAll,3));
           if ~isempty(BA)
               settings.channels = BA;
           else
               settings.channels = 1:size(settings.Spect.SpectAll,3);
           end
       end
   end
end

function Locs = load_locs(Locs,Spect)
   if ~isempty(Spect) & isfield(Spect,'SpectAll')
       numchans = size(Spect.SpectAll,3);
       Locs = lab_load_locs(Locs,[],numchans);
   else
       Locs = lab_load_locs(Locs);
   end
end

function settings = set_epochs(settings)
   if strcmp(settings.epochs,'single')
       settings.epochsnumber = 1;
   else
       settings.epochsnumber = [];
   end
end

function settings = set_channels(settings)
   if isempty(settings.Spect)
       return
   end
   if isfield(settings.Spect,'channels')
       settings.channels = lab_load_background(settings.channels,settings.Spect.channels,settings.Locs);
   else
       if ~isfield(settings.Spect,'names')
           if isfield(settings.Locs,'labels')
               plot.indexed = settings.channels;
               plot.Color = [1 1 1];
               plot.ColorIdx = [1 0 0];
               plot.LOCS = settings.Locs;
               plot.Title = 'Channels';
               settings.channels = lab_plot_locs(plot,1,0,0);
           else
               for j = 1:size(settings.Spect.SpectALL,3)
                   settings.Spect.names{j,1} = ['Spectrum_' num2str(j)];
               end
               settings.channels = listdlg('PromptString','Select Spectras','SelectionMode','multiple', ...
                   'Name','Select Spectras','ListString',settings.Spect.names,'InitialValue',settings.channels, ...
                   'CancelString','None','ListSize',[200 260]);
           end
       end
   end
end