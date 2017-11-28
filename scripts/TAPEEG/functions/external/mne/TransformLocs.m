for i = 1:size(data,1)
    header.locs.x(1,i) = header.chan_settings(1,i).loc(1,1);
    header.locs.y(1,i) = header.chan_settings(1,i).loc(2,1);
    header.locs.z(1,i) = header.chan_settings(1,i).loc(3,1);
end