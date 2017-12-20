function x=itab_filter_l100(in,Fs)

low=[104];

[b,a]=cheby2(24,60,low*2/Fs);

x=filtfilt(b,a,in);
