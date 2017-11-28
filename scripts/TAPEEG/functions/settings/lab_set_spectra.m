function [cfg,skipprocessing] = lab_set_spectra(cfg)

skipprocessing = 0;

if ~exist('cfg','var') | isempty(cfg) | ~isfield(cfg,'FFT') | ~isfield(cfg.FFT,'winsize')
    cfg.FFT.winsize = 4;
    cfg.FFT.hanpercent = 20;
    cfg.FFT.method = 'multitaper';
    cfg.FFT.lowfreq = 1;
    cfg.FFT.highfreq = 50;
    cfg.FFT.average = 'median';
    cfg.FFT.freqstep = 2;
end

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'FFT window (seconds)','winsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Hanning-window (percent)','hanpercent'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 100];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Method','method'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'welch','multitaper'};
Formats(end,1).callback = {@set_method,'@ALL','@ALL'};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Lowest frequency','lowfreq'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Highest frequency','highfreq'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Plot','average'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = {'mean','median'};

Prompt(end+1,:) = {'Interval frequency plot','freqstep'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 35;

[cfg.FFT,Cancelled] = inputsdlg(Prompt,'Spectral Analysis',Formats,cfg.FFT);
if isempty(cfg.FFT) | Cancelled == 1
    cfg.FFT = [];
    skipprocessing = 1;
    return
end

end

function settings = set_method(settings)
   if strcmp(settings.method,'multitaper')
       settings.hanpercent = 0;
   else
       settings.hanpercent = 20;
   end
end