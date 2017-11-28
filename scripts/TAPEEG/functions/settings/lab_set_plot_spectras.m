function [settings,skipprocessing] = lab_set_plot_spectras(settings)

disp ('   Ask for collect spectraldata settings')

skipprocessing = 0;

if ~exist('settings','var') | ~isfield(settings,'source')
    settings.file = '';
    settings.Spect = [];
    settings.mappingBA = [];
    settings.lowfreqpeak = 4;
    settings.highfreqpeak = 14;
    settings.lowfreqcog = 4;
    settings.highfreqcog = 14;
    settings.source = 'median';
    settings.MinPeak2Min = 1.3;
    settings.spectralbands = lab_get_spectralbands;
    settings.correctpf = true;
    settings.calcsingle = true;
    settings.qualityrange = 0.5;
    settings.qualityplot = true;
    settings.plotchans = false;
    settings.QUALITY = [];
    settings.mappings = [];
    settings.writeresults = false;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Spectra-File','file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.mat;*.xls;*.xlsx','Spectra-File (*.mat,*.xls,*.xlsx)'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 2];
Formats(end,1).callback = {@load_spectra,'@ALL','@ALL'};

Prompt(end+1,:) = {'Spectra','Spect'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_spectra,'@ALL','@ALL'};

Prompt(end+1,:) = {'Background Channels','mappingBA'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf];
Formats(end,1).enable = 'inactive';
Formats(end,1).callback = {@lab_load_background,'mappingBA','mappingBA',[],'Locs'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'LOCS','Locs'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_locs,'Locs','Locs','Spect'};

Prompt(end+1,:) = {'Mappings','mappings'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@load_mappings,'mappings','@ALL'};

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Source','source'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'mean','median','single'};
Formats(end,1).callback = {@set_source,'@ALL','@ALL'};

Prompt(end+1,:) = {'Spectral Bands','spectralbands'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).size = 70;
Formats(end,1).span = [1 2];
Formats(end,1).callback = {@lab_table_dialog,'spectralbands','spectralbands', ...
                             {'lowfreq','highfreq'},'Spectral Bands',1};

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Range peak frequency (Hz)    low:','lowfreqpeak'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'high:','highfreqpeak'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Range median frequency (Hz) low:','lowfreqcog'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;
Formats(end,1).span = [1 2];

Prompt(end+1,:) = {'high:','highfreqcog'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Minimal Peak2Min','MinPeak2Min'};
Formats(end+1,1).type = 'check';
Formats(end,1).callback = {@set_MinPeak2Min,'MinPeak2Min','MinPeak2Min'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Calculate single channel peaks','calcsingle'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];
Formats(end,1).format = 'integer';
Formats(end,1).callback = {@check_single,{'calcsingle','calcL2Rratio'},'calcsingle','calcL2Rratio'};

Prompt(end+1,:) = {'Manually correct peak frequency','correctpf'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Quality','qualityplot'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];

Prompt(end+1,:) = {'Range (+/-Hz)','qualityrange'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 30;

Prompt(end+1,:) = {'Evaluate','QUALITY'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).callback = {@get_quality,'QUALITY','@ALL'};

Prompt(end+1,:) = {' ',''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Write results','writeresults'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Plot spectras for channels','plotchans'};
Formats(end+1,1).type = 'check';
Formats(end,1).size = [-1 -1];
Formats(end,1).span = [1 3];

[settings,Cancelled] = inputsdlg(Prompt,'Plot spectra',Formats,settings);
if isempty(settings) | Cancelled == 1
    settings = [];
    skipprocessing = 1;
    return
end

end

function [calcsingle,calcL2Rratio] = check_single(calcsingle,calcL2Rratio)
  if calcsingle == false
      calcsingle = true;
  else
      calcsingle = false;
      calcL2Rratio = false;
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
               end
           case {'xls','xlsx'}
               try %#ok<TRYNC>
                   result = lab_read_xls(settings.file);
                   if ischar(result{1,1}) & (~isnumeric(result{1,2}) | ...
                           strcmp(result{1,1},'Patient') | strcmp(result{1,1},'Region') | ...
                           strcmp(result{1,1},'Channel'))
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
                               channels{j,1} = num2str(channels{j,1});
                           end
                       end
                       result = cell2mat(result(2:end,2:end));
                   elseif ischar(result{1,1}) & isnumeric(result{1,2})
                       channels = result(1:end,1);
                       result = cell2mat(result(1:end,2:end));
                       settings.lowfreq = 1;
                       settings.freqbin = 0.25;
                       Formats.format = 'float';
                       Formats.size = 40;
                       Formats(2,1).format = 'float';
                       Formats(2,1).size = 40;
                       settings = inputsdlg({'Lowest Frequency','lowfreq';'Freqbin','freqbin'}, ...
                           'Frequency range',Formats,settings,2);
                       freqlabel = settings.lowfreq:settings.freqbin:ceil(settings.freqbin*size(result,2));
                       freqlabel = freqlabel(1,1:size(result,2));
                   else
                       result = cell2mat(result);
                       for j = 1:size(result,1)
                           channels{j,1} = ['Spectrum_' num2str(j)];
                       end
                       settings.lowfreq = 1;
                       settings.freqbin = 0.25;
                       Formats.format = 'float';
                       Formats.size = 40;
                       Formats(2,1).format = 'float';
                       Formats(2,1).size = 40;
                       settings = inputsdlg({'Lowest Frequency','lowfreq';'Freqbin','freqbin'}, ...
                           'Frequency range',Formats,settings,2);
                       freqlabel = settings.lowfreq:settings.freqbin:ceil(settings.freqbin*size(result,2));
                       freqlabel = freqlabel(1,1:size(result,2));
                   end
                   settings.Spect.SpectAllF = freqlabel;
                   settings.Spect.SpectAllMean = result;
                   settings.Spect.SpectAllMedian = result;
                   settings.Spect.SpectAll = permute(result,[3 2 1]);
                   settings.Spect.channels = channels;
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
               settings.mappingBA = BA;
           else
               settings.mappingBA = 1:size(settings.Spect.SpectAll,3);
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

function Mappings = load_mappings(settings)
   if strcmp(settings.source,'single')
       Mappings = [];
       return
   end
   if isfield(settings.Spect,'SpectAll')
       settings.EXTRA.numdatachans = size(settings.Spect.SpectAll,3);
   elseif isfield(settings.Locs,'x')
       settings.EXTRA.numdatachans = size(settings.Locs.x,2);
   else
       settings = [];
   end
   Mappings = lab_load_mappings(settings.mappings,settings,[],settings.Locs);
end

function MinPeak2Min = set_MinPeak2Min(MinPeak2Min)
  if isempty(MinPeak2Min)
      MinPeak2Min = 1.3;
  end
  settings.MinPeak2Min = MinPeak2Min;
  
  Prompt = {'Minimal Peak2Min for detection of peak frequency','MinPeak2Min'};
  Formats.type = 'edit';
  Formats.format = 'float';
  Formats.limits = [0 inf];
  Formats.size = 30;
  
  [settings,Cancelled] = inputsdlg(Prompt,'Minimal Peak2Min',Formats,settings);
  if isempty(settings) | Cancelled == 1
      MinPeak2Min = 1.3;
  else
      MinPeak2Min = settings.MinPeak2Min;
  end
  
end

function settings = set_source(settings)
  if strcmp(settings.source,'single')
      settings.mappingBA = 1;
      settings.mappings = [];
  end
end

function QUALITY = get_quality(settings)
  doQUALITYchans = true;
  doQUALITYact = true;
  doQUALITYepoch = true;
  if isfield(settings.Spect,'header') & ~isfield(settings.Spect.header,'badchans')
      doQUALITYchans = false;
  end
  if isfield(settings.Spect,'header') & ~isfield(settings.Spect.header,'activationsexcluded')
      doQUALITYact = false;
  end
  if isfield(settings.Spect,'header') & ~isfield(settings.Spect.header,'quality')
      doQUALITYepoch = false;
  end
  QUALITY = lab_get_QUALITY(settings.QUALITY,0,doQUALITYchans,doQUALITYact,doQUALITYepoch);
end
