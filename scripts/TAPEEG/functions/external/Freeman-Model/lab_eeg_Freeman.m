% Calculates eeg signal with Freeman neural mass model
%
% [data,header] = lab_eeg_Freeman(chans,numtf,fs,matrix,lag)
%
% chans        = number of channels in result
% numtf        = number of timeframes in result
% fs           = sampling rate
% cmatrix      = matrix with connections (chans x chans)
% lag          = lag for connections in timeframes
% settings     = optional (see lab_set_eeg_Freeman)
%
% written by F. Hatz Vumc 2013

function [data,header,settings] = lab_eeg_Freeman(chans,numtf,fs,settings)

disp('    Calculate eeg data (Freeman neuronal mass model)')

if ~exist('settings','var')
    [settings,skipprocessing] = lab_set_eeg_Freeman([],1);
    if skipprocessing == 1
        settings = [];
        data = [];
        header = [];
        return
    else
        pause(0.2);
        chans = settings.chans;
        numtf = settings.numtf;
        fs = settings.fs;
    end
elseif isempty(settings)
    disp('   empty settings, set to standard values')
    settings.ae    = 25;       % rise time e s-1
    settings.be    = 175;      % decay time e s-1
    settings.ai    = 75;       % rise time i s-1
    settings.bi    = 225;      % decay time i s-1
    settings.vd    = 15;       % vthreshold mV
    settings.Qmax  = 250;      % max firing rate s-1
    settings.kei   = 5.0;      % coupling i to e population
    settings.kie   = 0.2;      % coupling e to i population
    settings.sigma = 1;        % stde
    settings.tref  = 5;        % refractory period
end
ae = settings.ae * 1000 / fs;
be = settings.be * 1000 / fs;
ai = settings.ai * 1000 / fs;
bi = settings.bi * 1000 / fs;
vd = settings.vd;
Qmax = settings.Qmax * 1000 / fs;
kei = settings.kei;
kie = settings.kie;
sigma = settings.sigma;
tref = settings.tref * 1000 / fs;
h     = 1/fs;  % time step
tf    = ceil((numtf+100)/tref)*tref;

fprintf('    ')
data = zeros(chans,tf);
for nchan = 1:chans
    ve    = zeros(1,2); % initial condition e cells
    dve   = zeros(1,2); % initial condition e cells, derivative
    Pe(1) = 0;          % initital condition Input e cells, sigmoid function
    vi    = zeros(1,2); % initial condition i cells
    dvi   = zeros(1,2); % initial condition i cells, derivative
    Pi(1) = 0;          % initial condition Input i cells, sigmoid function
    rng('default');
    rng('shuffle');
    % input noise
    tijdsmomenten = randi(tf-1,tf/tref,1);
    Iext = zeros(3,1);
    for i=2:tf
        if sum(i==tijdsmomenten)>0 && Iext(i-1,1)==0;
            Iext(i) = randn*20;
        else Iext(i) = 0;
        end
    end
    
    % freeman neural mass
    for i=2:tf
        % calculate membrane potential E population
        dve(i,:) = [ve(i-1,2) ae*be*Pe(i-1)-(ae+be)*ve(i-1,2)-(ae*be)*(ve(i-1,1))];
        ve(i,:) = ve(i-1,:) + dve(i,:)*h;
        % calculate membrane potential I population
        dvi(i,:) = [vi(i-1,2) ai*bi*Pi(i-1)-(ai+bi)*vi(i-1,2)-(ai*bi)*(vi(i-1,1))];
        vi(i,:) = vi(i-1,:) + dvi(i,:)*h;
        % potential to pulse with sigmoid function
        Pe(i) = kei*Qmax/(1+exp(-((vi(i,1) - vd)/sigma)))+Iext(i); %#ok<AGROW>
        Pi(i) = kie*Qmax/(1+exp(-((ve(i,1) - vd)/sigma))); %#ok<AGROW>
    end
    data(nchan,:) = ve(:,1);
    clearvars dve ve dvi Pe Pi
    fprintf('.')
end
disp(':')
data = data(:,101:numtf+100);

header.numtimeframes = size(data,2);
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.samplingrate = fs;
header.EEG_file = ['Freeman-model' num2str(size(data,1)) 'x' num2str(size(data,2))];
header.EEG_filepath = [];
header.channels = num2str((1:size(data,1))');
header.numauxchannels = 0;
tmp = clock;
header.year = tmp(1);
header.month = tmp(2);
header.day = tmp(3);
header.hour = tmp(4);
header.minute = tmp(5);
header.second = floor(tmp(6));
header.millisecond = (tmp(6) - floor(tmp(6)))*1000;
header.ref_chan = 'none';
