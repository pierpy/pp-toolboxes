function time_domain_features = calc_time_domain_features(epoch)

%Curve Length**************************************************************
depoch=epoch(2:end)-epoch(1:end-1);
% dy=0.04;%sample time, sample freq is 25Hz
curve_length=sum(abs(depoch));

%number of maxima and minima***********************************************
Number_of_min_and_max = length(findpeaks(epoch-mean(epoch),'THRESHOLD',4e-4));

%Root mean squared amplitude***********************************************
% rms = sqrt(mean(epoch.^2));

%Hjorth parameters*********************************************************
[ACTIVITY, MOBILITY, COMPLEXITY]=hjorth(epoch,0);

%zero crossings epoch, d en dd
zero_crossings = length(crossing(epoch-mean(epoch)));
dzero_crossings = length(crossing(depoch-mean(depoch)));
depoch=depoch-mean(depoch);
ddzero_crossings = length(crossing(depoch(2:end)-depoch(1:end-1)));

%Autoregressive modelling errors
Autoregressive_modelling_errors=zeros(1,9);
L=length(epoch);
for i=1:9
     a = armcov(epoch(1:L/2),i);
    est_epoch = filter([0 -a(2:end)],1,epoch(1+L/2:end));
    Autoregressive_modelling_errors(i) = sum(abs(est_epoch-epoch(1+L/2:end)));
end


%Skewness
Skewness = skewness(epoch);

%kurtosis
Kurtosis = kurtosis(epoch);

%Non-linear energy
Non_linear_energy = mean(NLEO(epoch)/length(epoch));

%variance epoch d en dd
variance = var(epoch);
dvariance = var(depoch);
ddvariance = var(depoch(2:end)-depoch(1:end-1));

%RMS amplituede
RMS_amplitude = sqrt((1/length(epoch)*sum(epoch.^2)));
%collect features**********************************************************

time_domain_features=[curve_length, Number_of_min_and_max, ...
    ACTIVITY, MOBILITY, COMPLEXITY, zero_crossings, dzero_crossings,...
    ddzero_crossings, Autoregressive_modelling_errors, Skewness,...
    Kurtosis, Non_linear_energy, variance, dvariance, ddvariance, RMS_amplitude];
