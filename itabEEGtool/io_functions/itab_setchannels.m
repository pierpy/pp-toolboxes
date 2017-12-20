function list_ch = itab_setchannels(header, code, quality)

%%%%%%%%%%% Find a list of channels (code)
%%%%%%%%%%% quality: 'working', 'unworking', 'all'

if nargin < 3
   quality = 'working';
end

switch code
    case 'ELEC'
        ch_type = 1;
    case 'MAG'
        ch_type = 2;
    case 'ELEC_REF'
        ch_type = 4;
    case 'MAG_REF'
        ch_type = 8;
    case 'AUX'
        ch_type = 16;
    case 'PARAM'
        ch_type = 32;
    case 'DIGIT'
        ch_type = 64;
end

nchan=header.nchan;
list_ch = [];
for i = 1:nchan
    switch quality
        case 'working'
            if (header.ch(i).type == ch_type) & (header.ch(i).flag == 0)
                list_ch = [list_ch i];
            end
        case 'unworking'
            if (header.ch(i).type == ch_type) & (header.ch(i).flag == 1)
                list_ch = [list_ch i];
            end
        case 'all'
            if header.ch(i).type == ch_type
                list_ch = [list_ch i];
            end
    end
end
   