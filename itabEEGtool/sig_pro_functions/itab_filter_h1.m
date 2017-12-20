function x=itab_filter_h1(in,Fs)

high=[0.7];

[b,a] = cheby2(5,40,high*2/Fs,'high');

x=filtfilt(b,a,in);