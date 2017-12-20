function x=itab_filter_h03(in,Fs)

% high=[0.09];
high=[0.07];

[b,a] = cheby2(4,40,high*2/Fs,'high');

x=filtfilt(b,a,in);
