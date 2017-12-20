function x=itab_filter_l80(in,Fs)

low=[84];

[b,a]=cheby2(22,60,low*2/Fs);

x=filtfilt(b,a,in);