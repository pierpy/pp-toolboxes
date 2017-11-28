function frequency_domain_features = calc_frequency_domain_features(epoch,EEG_power,f)

df=f(2)-f(1);

%Total Power 0-12Hz********************************************************
Total_Power=sum(EEG_power);

%Peak_frequency************************************************************
Peak_frequency=f(find(EEG_power==max(EEG_power((find(abs(f-10)==min(abs(f-10))):find(abs(f-32)==min(abs(f-32)))))))); %10==0.5 Hz, 526 == 32 Hz


%Spectral edge frequencies 80,90 en 95%
i=1;
while sum(EEG_power(1:i))<0.8*Total_Power
    i=i+1;
end
Spect_E_freq_80 = f(i);

i=1;
while sum(EEG_power(1:i))<0.9*Total_Power
    i=i+1;
end
Spect_E_freq_90 = f(i);

i=1;
while sum(EEG_power(1:i))<0.95*Total_Power
    i=i+1;
end
Spect_E_freq_95 = f(i);

%Power in 2 Hz subbands****************************************************
Power(1)=sum(EEG_power(2:floor(2/df)));
for i=2:11
    Power(i)=sum(EEG_power(floor((i-1)/df):floor((i+1)/df)));
end
%Normalized Power in 2 Hz subbands*****************************************

norm_Power=Power/Total_Power;

%Wavelet energy in 1-2Hz band => 5e coefficient****************************
%coefs = cwt(epoch,5,'db4');
%Wavelet_Energy = sum((1/25)*(coefs.^2));
Wavelet_Energy =0;
%collect features
frequency_domain_features=[Total_Power, Peak_frequency,Spect_E_freq_80,...
    Spect_E_freq_90, Spect_E_freq_95, Power,norm_Power, Wavelet_Energy];
