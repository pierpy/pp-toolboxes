function x = itab_filter_notch(in,notch,Fs)


if notch > 1000/Fs

    power = notch + notch*[-4 4]/Fs +[-0.3 0.3];
    [b,a] = cheby2(3,40,power*2/Fs,'stop');
    x = filter(b,a,in);

else
    x = in;
end