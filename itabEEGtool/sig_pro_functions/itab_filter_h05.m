function x=itab_filter_h05(in,Fs)

high=[0.35];

[b,a] = cheby2(4,40,high*2/Fs,'high');

x=filtfilt(b,a,in);