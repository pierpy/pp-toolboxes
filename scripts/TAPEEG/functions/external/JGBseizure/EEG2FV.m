function feature_vector = EEG2FV(EEGData,Fs)
%load('bandpass_05_32.mat')
%filtered_EEG=filtfilt(SOS,G,EEGData);
L=length(EEGData);%
w = hanning(L,'periodic');
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);
Y = fft(w.*EEGData,NFFT);
EEG_power=abs(Y(1:NFFT/2+1)).^2;
abs_volt = mean(abs(EEGData));
piek_piek = max(EEGData)-min(EEGData);
pre_proces_features = [abs_volt piek_piek];
Freq_features=calc_frequency_domain_features(EEGData,EEG_power,f);
time_domain_features = calc_time_domain_features(EEGData);
information_theory_features = calc_information_theory_features(EEGData,EEG_power,f);
speach_features=calc_speach_features(EEGData,w,Fs)';
feature_vector=[Freq_features time_domain_features...
    information_theory_features speach_features pre_proces_features];


       