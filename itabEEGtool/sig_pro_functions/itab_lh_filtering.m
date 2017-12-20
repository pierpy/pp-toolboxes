function sig_out = itab_lh_filtering(sig, h_filtering, l_filtering, Fs);

sig_out=sig;
switch h_filtering
    case 0.03
        sig_out = itab_filter_h03(sig,Fs);
    case 0.5
        sig_out = itab_filter_h05(sig,Fs);
    case 1
        sig_out = itab_filter_h1(sig,Fs);
    case 40
        sig_out = itab_filter_h40(sig,Fs);
end
sig = sig_out;
switch l_filtering
    case 80
        sig_out = itab_filter_l80(sig,Fs);
    case 100
        sig_out = itab_filter_l100(sig,Fs);
    case 150
        sig_out = itab_filter_l150(sig,Fs);
end
