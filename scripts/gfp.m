function [gfp] = globalFieldPower(data)
%global field power at time frame
% input: dati eeg [nch x 1]
% output: gfp --> global field power
%         pk  --> local maxima values
%         lk  --> local maxima indexs

u_hat = mean(data);
u = data - u_hat;
gfp=std(u);

