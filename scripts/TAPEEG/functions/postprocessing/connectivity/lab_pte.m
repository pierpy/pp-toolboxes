function [pte,settings] = lab_pte(data,header,settings)

pte = []; %#ok<NASGU>
if ~exist('settings','var')
    settings = [];
end

% filter data
[data,header] = lab_filtering(data,header);

% reduce to data channels
if isfield(header,'numdatachannels')
    [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels);
end

if ~isfield(settings,'phaseestimate')
    settings.phaseestimate = 'hilbert';
    settings.freqwindow = size(data,2);
    settings.usehann = false;
    settings.delta = 2;
    settings.deltafreq = true;
    settings.binfactor = 1;
end

if strcmp(settings.phaseestimate,'hilbert')
    phase = lab_hilbert(data,header,settings);
else
    disp('Abort: no valid method selected for PTE calculation');
end

Nchans = size(phase,1); % number of channels
Ntr = size(phase,2); % number of samples (time points)

if isfield(settings.delta,'deltafreq') & settings.deltafreq == false
    delta = settings.delta;
elseif isfield(header,'lowpass') & ~isempty(header.lowpass) & ...
        isfield(header,'highpass') & ~isempty(header.highpass)
    tmp = header.samplingrate / ((header.lowpass + header.highpass) / 2);
    delta = round(tmp / settings.delta);
    clearvars tmp
else
    % find optimal delta (based on zero crossings). delta: analysis lag
    if ~isfield(settings,'delta')
        settings.delta = 2;
    end
    tmp = [];
    for i = 1:size(phase,1)
        tmp = cat(2,tmp,diff(find(abs(diff(sign(phase(i,:))))==2)));
    end
    delta = round((2*mean(tmp)) / settings.delta);
    clearvars tmp
end
if delta < 1
    delta = 1;
end

% transform phase to ordinal numbers
H = zeros(Nchans,1);
for i = 1:Nchans
    % bin width according to Scott, 1992
    H(i,1) = (3.5 * circ_std(phase(i,:)')) / (Ntr^(1/3));
end
H = median(H);
H = round((2*pi/H)/2)*2; % H: number of bins
if isfield(settings,'binfactor') & ~isempty(settings.binfactor)
    H = round(H * settings.binfactor);
end
phase = round((phase + pi)/(2*pi) * (H - 1)) + 1; % phase in {1, 2, ..., H}

% set step
if ~isfield(settings,'step') | isempty(settings.step)
    settings.step = 0;
end
if settings.step > 0 & isfield(header,'samplingrate') & ~isempty(header.samplingrate)
    settings.step = settings.step * header.samplingrate;
else
    settings.step = size(phase,2) - delta;
end
if (settings.step + delta) > Ntr
    settings.step = size(phase,2) - delta;
end
Maxstep = floor((size(phase,2) - delta) / settings.step);

% calculate pte
pte = zeros(Nchans,Nchans,Maxstep);
for Nstep = 1:Maxstep
    phase1 = phase(:,1+delta+(Nstep-1)*settings.step:delta+Nstep*settings.step);
    phase2 = phase(:,1+(Nstep-1)*settings.step:Nstep*settings.step);
    L = size(phase2,2);
    phaseT = (phase1-1)*H + phase2;
    P1 = histc(phaseT,1:H^2,2) / L;
    R1 = - nansum(P1 .* log(P1),2);
    P3 = histc(phase2,1:H,2) / L;
    R3 = - nansum(P3 .* log(P3),2);
    R2 = zeros(Nchans,Nchans);
    R4 = zeros(Nchans,Nchans);
    for i = 1:Nchans
        P2 = histc((phase2-1)*H + repmat(phase2(i,:),Nchans,1),1:H^2,2) / L;
        R2(:,i) = - nansum(P2 .* log(P2),2);
        P4 = histc((phaseT-1)*H + repmat(phase2(i,:),Nchans,1),1:H^3,2) / L;
        R4(:,i) = - nansum(P4 .* log(P4),2);
    end
    pte(:,:,Nstep) = repmat(R1,[1 Nchans]) + R2 - repmat(R3,[1 Nchans]) - R4;
end

clearvars P1 P2 P3 P4 R1 R2 R3 R4 Nstep Maxstep i phase1 phase2 phaseT

for i = 1:size(pte,3)
    tmp = pte(:,:,i);
    tmp(1:Nchans+1:end) = 0;
    pte(:,:,i) = tmp;
end

end