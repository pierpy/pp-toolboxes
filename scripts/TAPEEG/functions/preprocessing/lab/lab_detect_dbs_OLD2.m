function [data,header,cfg,skipprocessing] = lab_detect_dbs_low(data,header,cfg)

disp('   detect stimulator artifacts (low voltage)')

skipprocessing = 0;

if ~exist('cfg','var') | ~isfield(cfg,'STIM') | ~isfield(cfg.STIM,'DETECT') |...
         ~isfield(cfg.STIM.DETECT,'stimrange') | isempty(cfg.STIM.DETECT.stimrange)
    cfg.STIM.DETECT.stimrange = [100,150;80,100;150,250];
end
if ~isfield(cfg.STIM.DETECT,'stimmaxpercent') | isempty(cfg.STIM.DETECT.stimmaxpercent)
    cfg.STIM.DETECT.stimmaxpercent = 70;
end
if ~exist('header','var')
    header= lab_create_header(data);
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = size(data,1);
end

if header.samplingrate < 500
    disp('     Abort: Samplingrate < 500Hz')
    skipprocessing = 1;
    return
end

fprintf('     do spectral analysis')
winsize = 2 * header.samplingrate;
hanning = hann(ceil(winsize * 0.4));
overlap = 0.2;
window = ones(1,winsize);
window(1:ceil(overlap*winsize)) = hanning(1:ceil(overlap*winsize));
window(end-ceil(overlap*winsize)+1:end) = hanning(end-ceil(overlap*winsize)+1:end);
window = window';
spectra = zeros(header.numdatachannels,floor(winsize/2)+1);
Counter = floor(header.numdatachannels/20);
for i = 1:header.numdatachannels
    [spectra(i,:),freqs] = pwelch(detrend(data(i,:)),window,1,winsize,header.samplingrate);
    if mod(i,Counter) == 0
        fprintf('.')
    end
end
disp(':')
freqs = freqs(1:floor(length(freqs) - length(freqs)/20));
spectra = spectra(:,1:length(freqs));
for i = 1:size(cfg.STIM.DETECT.stimrange,1)
    range(1,:) = [find(freqs>=cfg.STIM.DETECT.stimrange(i,1),1,'first') find(freqs>=cfg.STIM.DETECT.stimrange(i,2),1,'first')];
end
Freqs = zeros(size(range,1),2);
for i = 1:size(range,1)
    FreqT = zeros(1,size(spectra,1));
    AmpT = zeros(1,size(spectra,1));
    for j = 1:size(spectra,1)
        Idx = find(spectra(j,range(i,1):range(i,2)) == max(spectra(j,range(i,1):range(i,2))));
        Idx = Idx + range(1) - 1;
        FreqT(1,j) = freqs(Idx);
        AmpT(1,j) = max(spectra(j,range(i,1):range(i,2)));
    end
    Idx2 = find(AmpT > mean(AmpT)+std(AmpT));
    Freqs(i,1) = median(FreqT(Idx2));
    Freqs(i,2) = median(AmpT(Idx2));
    clearvars FreqT AmpT
end
if Freqs(1,2) >= max(Freqs(:,2)) * cfg.STIM.DETECT.stimmaxpercent/100
    Freq = Freqs(1,1);
else
    Freq = Freqs(Freqs(:,2)==max(Freqs(:,2)),1);
end
disp(['     set stimulator frequency to ' num2str(Freq) 'Hz'])

markerlength = floor(header.samplingrate/Freq);
MARK = round(1:header.samplingrate/Freq:size(data,2));
MARK = MARK(MARK+markerlength<=size(data,2));
events.POS = int64(MARK(:)');
events.DUR = int64(ones(1,length(MARK)));
events.OFF = int64(zeros(1,length(MARK)));
events.TYP = repmat({'STIM'},1,length(MARK));
header = lab_mix_markers(header,events);

end