function x=itab_filter_l40(in,Fs)

low=[43];

[b,a]=cheby2(12,40,low*2/Fs);

x=filtfilt(b,a,in);