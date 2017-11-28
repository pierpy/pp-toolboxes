function signal = lab_create_sinewave(freq,length,fs,randphase)

rng('default');
rng('shuffle');

if exist('randphase','var') & randphase == true
    phaselag = rand(1)*2*pi;
else
    phaselag = 0;
end
signal = sin((2*pi*freq/fs*(1:length*fs))+phaselag);