clear all
close all
clc


filename = 'D:\Pierpaolo\MATLAB\scripts\eeg fnirs\arnsim90_01\arnsim90_0101.edf';
parfilename = 'D:\Pierpaolo\MATLAB\scripts\eeg fnirs\arnsim90_01\arnsim90_0101.par';
locsfile = 'D:\Pierpaolo\MATLAB\scripts\microstatiCircadiano\OHBM\coord122.xyz';
performICA = 1;
saveEEGLABstr = 0;



[ EEG ] = itab_makeAnalysisStr( filename, parfilename, locsfile, saveEEGLABstr);
[ EEG ] = itab_skipintervals( EEG );
[ EEG ] = itab_runica( EEG );
[ EEG ] = itab_icaclassification( EEG );

sig = EEG.ICA.A*EEG.ICA.IC;
g_ic = EEG.ICA.gIC;
sig_clean = EEG.ICA.A(:,g_ic)*EEG.ICA.IC(g_ic,:);

for k=1:122
     plot(EEG.ICA.timeica, sig(k,:),'b', EEG.ICA.timeica,sig_clean(k,:),'r')
     title(num2str(EEG.chlabels{k}))
     pause
     close
end


[ EEG ] = itab_makeAnalysisStr( filename)
[ EEG ] = itab_makeAnalysisStr( filename, fulllocsfile, 0);
win = 1024*ceil(4/(ica_par.dec_fact+1));
spectral_estimator = spectrum.welch('Hamming',win,0);
for ix=1:length(gch)
    Power_Spectrum = psd(spectral_estimator,sig(ix,:),'NFFT',win,'Fs',Fsnew);
    mspettro(ix,:) = sqrt((Fsnew/win)*Power_Spectrum.Data)';
    Power_Spectrum = psd(spectral_estimator,sig_clean(ix,:),'NFFT',win,'Fs',Fsnew);
    mspettro_clean(ix,:) = sqrt((Fsnew/win)*Power_Spectrum.Data)';
    F = Power_Spectrum.Frequencies;
end

for k=1:length(gch)
     plot(F,mspettro(k,:),'b',F,mspettro_clean(k,:),'r')
     title(num2str(gch(k)))
     pause
end    
