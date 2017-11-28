% Calculates eeg signal with code from Green et al. 2012
%
% [data,header] = lab_eeg_ARmodel(chans,numtf,fs,coeff,matrix,lag)
%
% chans        = number of channels in result
% numtf        = number of timeframes in result
% fs           = samplingrate
% coeff        = coefficient to use (5 10 15 20 25 30)
% matrix       = matrix with connections (chans x chans)
% lag          = lag for connections in timeframes
%
% written by F. Hatz Vumc 2013
%
% original header:
% Hi Robert, I hope it results useful to you. The coefficients were obtained
% from the datasets in: 
%
% Halliday et al. "Using electroencephalography to study functional coupling
% between cortical activity and electromyograms during voluntary
% contractions in humans", Neuroscience Letters, 1998.
%
% Apart of the 10 coefficients used in the paper, I included sets of 5,
% 15, 20, 25 and 30 coefficients just for fun.

function [data,header] = lab_eeg_ARmodel(chans,numtf,fs,coeff)

rng('default');
rng('shuffle');

% set defaults
% arcoef5 = [1,-1.59897703764571,1.45489662416893,-1.18566923906890,0.712821225055818,-0.291111761067200];
arcoef10 = [1,-1.77579405983455,1.93696539859690,-2.00501691216994,1.84360254899545,-1.64161141314982,1.35349581410874,-1.05970608197994,0.754717657100847,-0.486341121350871,0.171993501544989];
% arcoef15 = [1,-1.84839630880226,2.15024944751414,-2.37787927881166,2.38522976934059,-2.34680529287637,2.21268249142518,-2.04604021480018,1.84491301896200,-1.62976023125658,1.31062298134843,-1.00089796292178,0.716183774795189,-0.451166018681652,0.252294736254125,-0.0902527845443400];
% arcoef20 = [1,-1.88310243058089,2.25094323705591,-2.55848415572674,2.66095687892344,-2.72599779377417,2.70277753882205,-2.65285110442062,2.56552195774275,-2.45591861515061,2.23087136989939,-1.99817303573475,1.76061287602921,-1.52338989333052,1.31534737145757,-1.10819409193818,0.889626193028041,-0.680873666322604,0.476633971911714,-0.295587322294684,0.120454214247162];
% arcoef25 = [1,-1.90052719590447,2.30414925062220,-2.65630398998198,2.80993519614930,-2.93020017537385,2.96325615425571,-2.96871831159711,2.93595965287671,-2.88142118879654,2.71274640209991,-2.53530964170299,2.35245432459760,-2.16367111682013,1.99321796502379,-1.81458635757152,1.61393055613000,-1.40885170383926,1.19583894446135,-0.985889359702215,0.756541564711037,-0.536356476978677,0.378351725935638,-0.229917177036260,0.127748697574389,-0.0434362321352398];
% arcoef30 = [1,-1.90052719590447,2.30414925062220,-2.65630398998198,2.80993519614930,-2.93020017537385,2.96325615425571,-2.96871831159711,2.93595965287671,-2.88142118879654,2.71274640209991,-2.53530964170299,2.35245432459760,-2.16367111682013,1.99321796502379,-1.81458635757152,1.61393055613000,-1.40885170383926,1.19583894446135,-0.985889359702215,0.756541564711037,-0.536356476978677,0.378351725935638,-0.229917177036260,0.127748697574389,-0.0434362321352398];
% error5 = 0.140711391147563;
error10 = 0.119412750998595;
% error15 = 0.112009394211359;
% error20 = 0.107678216570111;
% error25 = 0.106145289651938;
% error30 = 0.106145289651938;

% set selection to coeff
if ~exist('coeff','var')
    arcoef=arcoef10;
    generror=error10;
else
    eval(['arcoef=arcoef' num2str(coeff) ';'])
    eval(['generror=error' num2str(coeff) ';'])
end
if fs ~= 500
    numtf2 = numtf * 500 / fs;
else
    numtf2 = numtf;
end

if ~exist('numtf','var')
    numtf2 = 5000;
end
if ~exist('chans','var')
    chans = 1;
end

elimsam=3000; %Samples to eliminate after modelling
vecsize=round(numtf2); %Size of the desired segments
data = zeros(chans,vecsize);
for i = 1:chans    
    %Generation error
    varerror=sqrt(generror);
    %Generation of independent source
    
    ran1=randn(1,vecsize+elimsam);
    %Synthesizing, EEG
    aux=filter(1,arcoef,varerror*ran1,[],2);
    %Eliminating the first elisam samples
    data(i,:)=aux(elimsam+1:end);
end

if fs ~= 500
   data = lab_resample_data(data,500,fs);
end

header.numtimeframes = size(data,2);
header.numchannels = size(data,1);
header.numdatachannels = header.numchannels;
header.samplingrate = fs;
header.EEG_file = ['AR-model' num2str(size(data,1)) 'x' num2str(size(data,2))];
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

return
