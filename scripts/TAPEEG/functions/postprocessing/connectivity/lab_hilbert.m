% Helper script for lab_calculate_connect_phase
%   Calculate instantaneous phase for all channels using hilbert
%   transformation
%
% written by F.Hatz 2014

function [phase,settings] = lab_hilbert(data,header,settings,cfgfilt)

if ~exist('cfgfilt','var')
    cfgfilt = [];
end
if ~exist('settings','var') | ~isfield(settings,'freqwindow')
    settings.freqwindow = size(data,2);
end
if ~isfield(settings,'usehann')
    settings.usehann = false;
end

step = floor(settings.freqwindow/2);
stepnr = floor(size(data,2) / step);
stepmax = floor(stepnr - (settings.freqwindow/step) + 1);
if stepmax > 1
    datstart = floor(settings.freqwindow / 2) - ceil(step/2) + 1;
    datend = datstart + step -1 ;
else
    datstart = 1;
    datend = settings.freqwindow;
end

if stepmax < 1
    stepmax = 1;
end

if exist('cfgfilt','var') & stepmax > 1
    data = lab_filter(data,header,cfgfilt,'novrb');
    cfgfilt = [];
end

if settings.usehann == 1
    hanningwin = ones(1,settings.freqwindow);
    hanntmp = hann((ceil(settings.freqwindow/2)))';
    hannS = ceil(size(hanntmp,2) / 2);
    hanningwin(1,1:hannS) = hanntmp(1,1:hannS);
    hanningwin(1,end-hannS:end) = hanntmp(1,end-hannS:end);
end
clearvars hannS hanntmp



phase = [];
disp(['     calculate hilbert with ' num2str(stepmax) ' steps (window = ' num2str(settings.freqwindow) ' tf)'])
for nstep = 1:stepmax
    m = (nstep-1)*step + 1;
    n = (nstep-1)*step + settings.freqwindow;
    datatmp = data(:,m:n);
    
    % apply hanning window
    if exist('hanningwin','var')
        datatmp = datatmp .* repmat(hanningwin,size(data,1),1);
    end
 
    % do filtering if not already done
    if ~isempty(cfgfilt)
        cfgfilt.nodisp = 1;
        datatmp = lab_filter(datatmp,header,cfgfilt,'novrb');
    end
    
    % do hilbert and angle
    datahilbtmp = hilbert(datatmp')';
    if nstep == 1
        dataphase = angle(datahilbtmp(:,1:datend));
    elseif nstep == stepmax
        dataphase = angle(datahilbtmp(:,datstart:end));
    else
        dataphase = angle(datahilbtmp(:,datstart:datend));
    end
    phase = cat(2,phase,dataphase);
end
clearvars dataphase datahilbtmp datatmp m n