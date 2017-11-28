function [peackData, gfp] = globalFieldPower(data)
% input: dati eeg [nch x nsamples]
% output: gfp --> global field power
%         peackData --> data values at the gfp maxima

 
 for k=1:size(data,2)
     gfp(k)=std(data(:,k));
 end
 
 
 [pk,lk] = findpeaks(gfp);
 peackData = data(:,lk);
 
%  
% plot(gfp)
% hold on
% plot(lk,pk,'or')