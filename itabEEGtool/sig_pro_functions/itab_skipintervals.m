function [ EEG ] = itab_skipintervals( EEG )
%ITAB_SKIPINTERVALS Summary of this function goes here
%   Detailed explanation goes here
    newfs = EEG.srate;
    skipped_intervals = EEG.skipped;
    ntpdata = EEG.pnts;
    indata = EEG.data;
    offset=10*newfs;

    %------> set TIME AXIS <------%
    timeline_dec = 1/newfs*[1:round(ntpdata)];
    time_position = ones(size(timeline_dec));
    for z=1:size(skipped_intervals,1)
        time_position(timeline_dec >= skipped_intervals(z,1) & timeline_dec <= skipped_intervals(z,2)) = 0;
    end
    time_position(1: offset) = 0;
    time_position(end-offset+1: end) = 0;
    time_flag = find(time_position == 1); 
    time_ica = timeline_dec(time_flag); 
    dleng = length(indata);
    outdata = zeros(size(indata,1), dleng);
    outdata(:, time_flag) = indata(:, time_flag);
    icadata = indata(:, time_flag);
    EEG.ICA.icadata = icadata;
    EEG.data = outdata;
    EEG.timeflag = time_flag;
    EEG.time = timeline_dec;
    EEG.ICA.timeica = time_ica;
end

