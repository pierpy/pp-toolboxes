function [h_filtering, l_filtering, freq_notch, dec, ica_method, ica_iter, num_segments,...
    skipped_intervals, study, file_evt, ncond, event, low_freq, up_freq, ech, res_chan] = itab_readparfile(filename_par)

fp=fopen(filename_par,'rt');
line=fgetl(fp);
line=fgetl(fp);

stringa = textscan(fp, '%s %s %s %s %s %s %f', 2);
h_filtering = stringa{1,end}(1);
l_filtering = stringa{1,end}(2);

stringa = textscan(fp, '%s %s %s %s %s %s %s', 1);

freq_notch = [];
if(strcmp(stringa{1,end},'none'))
    freq_notch = [];
    stringa = textscan(fp, '%s', 1);
else
    freq_notch = str2double(stringa{1,end});
    stringa = textscan(fp, '%s', 1);
    while(~strcmp(stringa{1,1},'Decimation'))
        freq_notch = [freq_notch str2double(stringa{1,end})];
        stringa = textscan(fp, '%s', 1);
    end
end

stringa = textscan(fp, '%s %f', 1);
dec = stringa{1,end};

stringa = textscan(fp, '%s %s %s', 1);
ica_method = cell2mat(stringa{1,end});

stringa = textscan(fp, '%s %s %s %f', 1);
ica_iter = stringa{1,end};

stringa = textscan(fp, '%s %s %s %s %s %s %f', 1);
num_segments = stringa{1,end};

skipped_intervals = zeros(num_segments,2);

if(num_segments~=0)
    stringa = textscan(fp, '%s %s %s',1); 
    for k = 1:num_segments
        skipped_intervals(k,1) = cell2mat(textscan(fp, '%f',1));
        skipped_intervals(k,2) = cell2mat(textscan(fp, '%f',1));
    end
% else
%    stringa = textscan(fp, '%s %s %s',1); 
end

stringa = textscan(fp, '%s %s %s',1);
tmp_str = stringa{1,end};
study = 0;
if(strcmp(tmp_str,'ongoing'))
        study = 0;
elseif(strcmp(tmp_str,'event-related'))
        study = 1;
elseif(strcmp(tmp_str,'block-design'))
        study = 2;
elseif(strcmp(tmp_str,'event-related-erd_ers'))
        study = 3;
end


ncond = 0;
file_evt = [];
event = [];
low_freq = [];
up_freq = [];

if study > 0
    stringa = textscan(fp, '%s %s %s %s',1);
    file_evt = stringa{1,end};
    file_evt =  file_evt{1,1};
    
    stringa = textscan(fp, '%s %s %s %d',1);
    ncond = stringa{1,end};
    
    for i=1:ncond
        stringa = textscan(fp, '%s %s %s',1);
        tmp_str = stringa{1,end}; %condition name
        event(i).name = tmp_str{1,1};
        
        if (study==1)|(study==3)
            stringa = textscan(fp, '%s %s %d',1);
            event(i).labels = stringa{1,end}; %trigger code
        else
            stringa = textscan(fp, '%s %s %s %d',1);
            event(i).labels_on = stringa{1,end};
            stringa = textscan(fp, '%s %s %s %d',1);
            event(i).labels_off = stringa{1,end};
        end
      
        stringa = textscan(fp, '%s %s %f',1);
        event(i).t_pre = stringa{1,end};
        stringa = textscan(fp, '%s %s %f',1);
        event(i).t_post = stringa{1,end};
        stringa = textscan(fp, '%s %s %s %f',1);
        event(i).t_baseline_on = stringa{1,end};
        stringa = textscan(fp, '%s %s %s %f',1);
        event(i).t_baseline_off = stringa{1,end};
    end

    if (study==2)|(study==3)
        stringa = textscan(fp, '%s %s %s %f',1);
        low_freq = stringa{1,end};
        stringa = textscan(fp, '%s %s %s %f',1);
        up_freq = stringa{1,end};
    end
    
end

% Read electrich channel list
stringa = textscan(fp, '%s %s %s',1)
ech = fscanf(fp, '%d')';

stringa = textscan(fp, '%s %s %s %s',1)
res_chan = stringa{1,end};
res_chan = res_chan{1,1};

fclose(fp)



