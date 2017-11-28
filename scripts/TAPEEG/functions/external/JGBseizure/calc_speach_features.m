function speach_features = calc_speach_features(x,w,fs)
x = x(:)';          % make row vector
w = w(:)';
%  x=x-mean(x);
nfft = 2^nextpow2(length(x));

 X = fft(x.*w,nfft);  
f = fs/2*linspace(0,1,nfft/2+1);
% plot(2*abs(X(1:nfft/2+1)),'r') 
n2 = 1+floor(nfft/2);
mfbk15 = melbankm(150,nfft,fs)/2;
mfbk15 = mfbk15(1:15,:);
% mfbk17 = melbankm(17,nfft,fs)/2;
mfbk20 = melbankm(200,nfft,fs)/2;
mfbk20 = mfbk20(1:20,:);

logFBE=log10(mfbk15*abs(X(1:n2))');%normalisatie! 

% relFBE=linFBE/sum(linFBE);%Relative filterbank energies


% FF=zeros(1,15);%Frequency-filtered band energies
% RSD=zeros(1,15);%Relative spectral differences

padded_logFBE=[0; log10(mfbk15*abs(X(1:n2))'); 0];
FF=padded_logFBE(3:end)-padded_logFBE(1:end-2);
%  FF=ff(2:end-1);
%  b=[-1 0 1];
%  FF=filter(b,1,padded_logFBE);
% padded_linFBE=[0 linFBE 0];
% for i=2:15+1
%    FF(i-1)= padded_logFBE(i)-(padded_logFBE(i+1)+padded_logFBE(i-1));
%    RSD(i-1)= (padded_linFBE(i+1)-padded_linFBE(i-1))/(padded_linFBE(i-1)+padded_linFBE(i)+padded_linFBE(i+1));
% end


%  cc = dct(log10(mfbk20*abs(X(1:n2))')/norm(log10(mfbk20*abs(X(1:n2))')));
%  cc=cc/norm(cc);
% plot(cc(2:end));hold on

  cc = dct(log10(mfbk20*(abs(X(1:n2)))')); 
%  plot(cc(2:end),'r')
%  cc = dct(log10(mfbk20*abs(X(1:n2))').^2); 
%  plot(cc(2:end),'g')
%cc = abs(fft(Xmel));
% cc=sort(abs(cc),'descend');
CC = cc(1:15);
% dCC = CC(2:end)-CC(1:end-1);
% CC_energies = idct(cc);
speach_features=[logFBE; FF; CC];
% speach_features=[linFBE relFBE logFBE FF RSD CC];