function x=itab_filter_h40(in,Fs)

high=[37];

[b,a] = cheby2(12,40,high*2/Fs,'high');

x=filtfilt(b,a,in);