function [pk,lk] = gfpPeack(gfp)
% input : global field power [1 x nframes]
[pk,lk] = findpeaks(gfp);



% eventuale plot
% plot(gfp)
% hold on
% plot(lk,pk,'or')
% EEG.chanlocs = readlocs('19ch2.xyz');
% topoplot(b(1,:), EEG.chanlocs)