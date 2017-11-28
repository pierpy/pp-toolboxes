% [data,header,cfg] = lab_read_fiff(filename,cfg) - read Elekta *.fiff
%
% In case of Elekta 306 channel file, channel will be reordered to
% 1. transversal gradiometers, 2. radial gradiometers 3. magnetometers
%
% In case filename contents 'part1of2_VE.fif', the second file
% 'part2of2_VE.fif' is added
%
% In case of 116 channels (AAL-atlas), you can optionally reduce channels
% to 78 regions (standard VUMC Amsterdam)
%
% written by F. Hatz 2012

function [data,header,cfg] = lab_read_fiff(filename,cfg,nodata,segment)

if ~exist('segment','var')
    segment = [];
end
if ~exist('nodata','var')
    nodata = false;
end
if ~exist('cfg','var')
    cfg = [];
end
if ~exist('filename','var')
    [filename,filepath] = uigetfile('*.fif;*.fiff','Select file');
    filename = fullfile(filepath,filename);
    clearvars filepath
end

warning off %#ok<WNOFF>
[hdr] = fiff_setup_read_raw(filename);
if nodata == false
    if ~isempty(segment) & length(segment) > 1
        segment = segment + hdr.firstsample - 1;
        if segment(1) > hdr.last_samp-1
            segment(1) = hdr.last_samp-1;
        end
        if segment(2) > hdr.last_samp
            segment(2) = hdr.last_samp;
        end
        [data,~] = fiff_read_raw_segment(hdr,segment(1),segment(2));
    else
        [data,~] = fiff_read_raw_segment(hdr);
    end
else
    data = zeros(hdr.info.nchan,1);
end
warning on %#ok<WNON>

header.year=0;
header.month=0;
header.day=0;
header.hour=0;
header.minute=0;
header.second=0;
header.millisecond=0;
header.numauxchannels = 0;
header.numchannels = size(data,1);
header.numtimeframes = size(data,2);
header.highpass = hdr.info.highpass;
header.lowpass = hdr.info.lowpass;
header.channels = char(hdr.info.ch_names');
header.chan_settings = hdr.info.chs;
header.samplingrate = hdr.info.sfreq;
header.begin = double(hdr.first_samp) / hdr.info.sfreq;
header.end = double(hdr.last_samp) / hdr.info.sfreq;
header.orig = hdr;

% create locs
for i = 1:size(data,1)
    header.locs.x(1,i) = header.chan_settings(1,i).loc(1,1)*1000;
    header.locs.y(1,i) = header.chan_settings(1,i).loc(2,1)*1000;
    header.locs.z(1,i) = header.chan_settings(1,i).loc(3,1)*1000;
    header.locs.labels{1,i} = header.chan_settings(1,i).ch_name;
end

try  %#ok<TRYNC>
    warning off %#ok<WNOFF>
    [grad] = mne2grad(hdr.info);
    warning on %#ok<WNON>
    if ~isempty(grad) & ~isempty(grad.coilpos)
        header.locs.grad = ft_convert_units(grad,'mm');
    end
end

header.datatype = 'meg';
header.ref_chan = 'none';

% Define data channels
header.numauxchannels = 0;
for i = 1:header.numchannels
    if strcmp(header.channels(i,1:3),'MEG')
        header.numdatachannels = i;
    elseif strcmp(header.channels(i,1:3),'STI') & ~isfield(header,'numdatachannels')
        header.numdatachannels = i - 1;
    end
end
if isfield(header,'numdatachannels')
    header.numauxchannels = header.numchannels - header.numdatachannels;
else
    header.numdatachannels = header.numchannels;
    header.numauxchannels = 0;
end

% Read events
header.events.POS = [];
header.events.DUR = [];
header.events.OFF = [];
header.events.TYP = [];
markerchan = [];
j = 1;
for i = header.numdatachannels+1:header.numchannels
    if strcmp(header.channels(i,1:4),'STI0')
        markerchan = [markerchan i]; %#ok<AGROW>
        markers{j,1} = header.channels(i,1:6); %#ok<AGROW>
        markers{j,2} = find(data(i,:) > 0); %#ok<AGROW>
        j = j + 1;
    elseif  strcmp(header.channels(i,1:3),'STI')
        markerchan = [markerchan i]; %#ok<AGROW>
    end
end
clearvars j
if exist('markers','var')
    for i = 1:size(markers,1)
        tmp = markers{i,2};
        if ~isempty(tmp)
            n = 1;
            tmp2 = tmp(1,1);
            tmp2(1,2) = 0;
            for j = 2:size(tmp,2)
                if tmp(1,j) == (tmp(1,j-1) + 1);
                    tmp2(n,2) = tmp2(n,2) + 1;
                else
                    n = n + 1;
                    tmp2(n,1) = tmp(1,j);
                    tmp2(n,2) = 0;
                end
            end
            tmp2(tmp2(:,2) == 0,2) = 1;
            for j = 1:size(tmp2,1)
                header.events.POS = [header.events.POS int64(tmp2(j,1))];
                header.events.DUR = [header.events.DUR int64(tmp2(j,2))];
                header.events.OFF = [header.events.OFF int64(0)];
                header.events.TYP = [header.events.TYP cellstr(markers{i,1})];
            end
            clearvars tmp tmp2 n j
        end
    end
    clearvars i markers
end
% Resort events
if size(header.events.POS,2) > 0
    resort = 1:size(header.events.POS,2);
    resort(2,:) = header.events.POS;
    resort = sortrows(resort',2)';
    resort = resort(1,:);
    header.events.POS = header.events.POS(1,resort);
    header.events.DUR = header.events.DUR(1,resort);
    header.events.OFF = header.events.OFF(1,resort);
    header.events.TYP = header.events.TYP(1,resort);
    clearvars resort
end
if isempty(header.events.POS)
    header = rmfield(header,'events');
end

% Eliminate marker-channels in data
if ~isempty(markerchan)
    includechannels = setdiff(1:header.numchannels,markerchan);
    [data,header] = lab_reduce_channels(data,header,includechannels);
end
clearvars markerchan tmp

% Locs to spherical coordinates
header.locs = lab_locs2sph(header.locs);

% Correct data
if strcmp(header.channels(1,1:3),'MEG') & header.numdatachannels == 306
    % Sort data: 1. transversal gradiometers, 2. radial gradiometers 3. magnetometers
    sensororder = {'MEG0113';'MEG0122';'MEG0132';'MEG0143';'MEG0213';'MEG0222'; ...
        'MEG0232';'MEG0243';'MEG0313';'MEG0322';'MEG0333';'MEG0343';'MEG0413'; ...
        'MEG0422';'MEG0432';'MEG0443';'MEG0513';'MEG0523';'MEG0532';'MEG0542'; ...
        'MEG0622';'MEG0633';'MEG0642';'MEG0713';'MEG0723';'MEG0733';'MEG0743'; ...
        'MEG0813';'MEG0822';'MEG0913';'MEG0923';'MEG0932';'MEG0942';'MEG1013'; ...
        'MEG1023';'MEG1032';'MEG1043';'MEG1112';'MEG1123';'MEG1133';'MEG1142'; ...
        'MEG1213';'MEG1223';'MEG1232';'MEG1243';'MEG1312';'MEG1323';'MEG1333'; ...
        'MEG1342';'MEG1412';'MEG1423';'MEG1433';'MEG1442';'MEG1512';'MEG1522'; ...
        'MEG1533';'MEG1543';'MEG1613';'MEG1622';'MEG1632';'MEG1643';'MEG1713'; ...
        'MEG1722';'MEG1732';'MEG1743';'MEG1813';'MEG1822';'MEG1832';'MEG1843'; ...
        'MEG1912';'MEG1923';'MEG1932';'MEG1943';'MEG2013';'MEG2023';'MEG2032'; ...
        'MEG2042';'MEG2113';'MEG2122';'MEG2133';'MEG2143';'MEG2212';'MEG2223'; ...
        'MEG2232';'MEG2242';'MEG2312';'MEG2323';'MEG2332';'MEG2343';'MEG2412'; ...
        'MEG2423';'MEG2433';'MEG2442';'MEG2512';'MEG2522';'MEG2533';'MEG2543'; ...
        'MEG2612';'MEG2623';'MEG2633';'MEG2642';'MEG0112';'MEG0123';'MEG0133'; ...
        'MEG0142';'MEG0212';'MEG0223';'MEG0233';'MEG0242';'MEG0312';'MEG0323'; ...
        'MEG0332';'MEG0342';'MEG0412';'MEG0423';'MEG0433';'MEG0442';'MEG0512'; ...
        'MEG0522';'MEG0533';'MEG0543';'MEG0612';'MEG0613';'MEG0623';'MEG0632'; ...
        'MEG0643';'MEG0712';'MEG0722';'MEG0732';'MEG0742';'MEG0812';'MEG0823'; ...
        'MEG0912';'MEG0922';'MEG0933';'MEG0943';'MEG1012';'MEG1022';'MEG1033'; ...
        'MEG1042';'MEG1113';'MEG1122';'MEG1132';'MEG1143';'MEG1212';'MEG1222'; ...
        'MEG1233';'MEG1242';'MEG1313';'MEG1322';'MEG1332';'MEG1343';'MEG1413'; ...
        'MEG1422';'MEG1432';'MEG1443';'MEG1513';'MEG1523';'MEG1532';'MEG1542'; ...
        'MEG1612';'MEG1623';'MEG1633';'MEG1642';'MEG1712';'MEG1723';'MEG1733'; ...
        'MEG1742';'MEG1812';'MEG1823';'MEG1833';'MEG1842';'MEG1913';'MEG1922'; ...
        'MEG1933';'MEG1942';'MEG2012';'MEG2022';'MEG2033';'MEG2043';'MEG2112'; ...
        'MEG2123';'MEG2132';'MEG2142';'MEG2213';'MEG2222';'MEG2233';'MEG2243'; ...
        'MEG2313';'MEG2322';'MEG2333';'MEG2342';'MEG2413';'MEG2422';'MEG2432'; ...
        'MEG2443';'MEG2513';'MEG2523';'MEG2532';'MEG2542';'MEG2613';'MEG2622'; ...
        'MEG2632';'MEG2643';'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211'; ...
        'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341'; ...
        'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531'; ...
        'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721'; ...
        'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931'; ...
        'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121'; ...
        'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311'; ...
        'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441'; ...
        'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631'; ...
        'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821'; ...
        'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011'; ...
        'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141'; ...
        'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331'; ...
        'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521'; ...
        'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641';};
    resort = zeros(1,header.numdatachannels);
    for i = 1:header.numdatachannels
        tmp = strfind(sensororder,header.locs.labels{1,i});
        resort(1,i) = find(~cellfun('isempty',tmp));
        clearvars tmp
    end
    if min(resort) > 0
        resort2 = [resort (header.numdatachannels+1:header.numchannels)];
        data(resort2,:) = data;
        header.chan_settings(1,resort2) = header.chan_settings(1,1:header.numchannels);
        header.channels(resort2,:) = header.channels;
        header.locs.x = header.locs.x(1,1:header.numdatachannels);
        header.locs.y = header.locs.y(1,1:header.numdatachannels);
        header.locs.z = header.locs.z(1,1:header.numdatachannels);
        header.locs.radius = header.locs.radius(1,1:header.numdatachannels);
        header.locs.theta = header.locs.theta(1,1:header.numdatachannels);
        header.locs.labels = header.locs.labels(1,1:header.numdatachannels);
        header.locs.x(1,resort) = header.locs.x;
        header.locs.y(1,resort) = header.locs.y;
        header.locs.z(1,resort) = header.locs.z;
        header.locs.radius(1,resort) = header.locs.radius;
        header.locs.theta(1,resort) = header.locs.theta;
        header.locs.labels(1,resort) = header.locs.labels;
        if isfield(header.locs,'grad')
            header.locs.grad.label(resort,1) = header.locs.grad.label;
            header.locs.grad.tra(resort,:) = header.locs.grad.tra;
        end
        clearvars resort resort2
        header.splitchans = [1,102;103,204;205,306];
    end
end

if strcmp(filename(end-14:end),'part1of2_VE.fif') & size(data,1) == 100
    % In case of splitted beamformer output add second file
    filename2 = [filename(1:end-15) 'part2of2_VE.fif'];
    if exist(filename2,'file')
        disp(['     add ' filename(1:end-15) 'part2of2_VE.fif'])
        cfg.filename = [filename(1:end-16) '.fif'];
        [hdr2] = fiff_setup_read_raw(filename2);
        [data2,~] = fiff_read_raw_segment(hdr2);
        if size(data2,1) == 17
            data = cat(1,data(1:end,:),data2(1:end-1,:));
            header.numchannels = size(data,1);
            header.orig.info.chs = [header.orig.info.chs(1,1:100) hdr2.info.chs];
            header.orig.info.ch_names = [header.orig.info.ch_names(1,1:100) hdr2.info.ch_names];
            header.chan_settings = [header.chan_settings(1,1:100) hdr2.info.chs];
            header.orig.info.nchan = 117;
            for i = 1:size(data,1)
                header.locs.x(1,i) = header.chan_settings(1,i).loc(1,1)*1000;
                header.locs.y(1,i) = header.chan_settings(1,i).loc(2,1)*1000;
                header.locs.z(1,i) = header.chan_settings(1,i).loc(3,1)*1000;
                header.locs.labels{1,i} = header.chan_settings(1,i).ch_name;
            end
            header.locs = lab_locs2sph(header.locs);
            Labels = {'Precentral_L';'Precentral_R';'Frontal_Sup_L';'Frontal_Sup_R'; ...
                'Frontal_Sup_Orb_L';'Frontal_Sup_Orb_R';'Frontal_Mid_L';'Frontal_Mid_R'; ...
                'Frontal_Mid_Orb_L';'Frontal_Mid_Orb_R';'Frontal_Inf_Oper_L'; ...
                'Frontal_Inf_Oper_R';'Frontal_Inf_Tri_L';'Frontal_Inf_Tri_R'; ...
                'Frontal_Inf_Orb_L';'Frontal_Inf_Orb_R';'Rolandic_Oper_L'; ...
                'Rolandic_Oper_R';'Supp_Motor_Area_L';'Supp_Motor_Area_R'; ...
                'Olfactory_L';'Olfactory_R';'Frontal_Sup_Medial_L';'Frontal_Sup_Medial_R'; ...
                'Frontal_Mid_Orb_L';'Frontal_Mid_Orb_R';'Rectus_L';'Rectus_R';'Insula_L'; ...
                'Insula_R';'Cingulum_Ant_L';'Cingulum_Ant_R';'Cingulum_Mid_L'; ...
                'Cingulum_Mid_R';'Cingulum_Post_L';'Cingulum_Post_R';'Hippocampus_L'; ...
                'Hippocampus_R';'ParaHippocampal_L';'ParaHippocampal_R';'Amygdala_L'; ...
                'Amygdala_R';'Calcarine_L';'Calcarine_R';'Cuneus_L';'Cuneus_R';'Lingual_L'; ...
                'Lingual_R';'Occipital_Sup_L';'Occipital_Sup_R';'Occipital_Mid_L'; ...
                'Occipital_Mid_R';'Occipital_Inf_L';'Occipital_Inf_R';'Fusiform_L'; ...
                'Fusiform_R';'Postcentral_L';'Postcentral_R';'Parietal_Sup_L'; ...
                'Parietal_Sup_R';'Parietal_Inf_L';'Parietal_Inf_R';'SupraMarginal_L'; ...
                'SupraMarginal_R';'Angular_L';'Angular_R';'Precuneus_L';'Precuneus_R'; ...
                'Paracentral_Lobule_L';'Paracentral_Lobule_R';'Caudate_L';'Caudate_R'; ...
                'Putamen_L';'Putamen_R';'Pallidum_L';'Pallidum_R';'Thalamus_L';'Thalamus_R'; ...
                'Heschl_L';'Heschl_R';'Temporal_Sup_L';'Temporal_Sup_R';'Temporal_Pole_Sup_L'; ...
                'Temporal_Pole_Sup_R';'Temporal_Mid_L';'Temporal_Mid_R';'Temporal_Pole_Mid_L'; ...
                'Temporal_Pole_Mid_R';'Temporal_Inf_L';'Temporal_Inf_R';'Cerebelum_Crus1_L'; ...
                'Cerebelum_Crus1_R';'Cerebelum_Crus2_L';'Cerebelum_Crus2_R';'Cerebelum_3_L'; ...
                'Cerebelum_3_R';'Cerebelum_4_5_L';'Cerebelum_4_5_R';'Cerebelum_6_L'; ...
                'Cerebelum_6_R';'Cerebelum_7b_L';'Cerebelum_7b_R';'Cerebelum_8_L'; ...
                'Cerebelum_8_R';'Cerebelum_9_L';'Cerebelum_9_R';'Cerebelum_10_L'; ...
                'Cerebelum_10_R';'Vermis_1_2';'Vermis_3';'Vermis_4_5';'Vermis_6';'Vermis_7'; ...
                'Vermis_8';'Vermis_9';'Vermis_10'};
            Label = {'LPrecent';'RPrecent';'LFroSup';'RFroSup';'LFroSupO';'RFroSupO'; ...
                'LFroMid';'RFroMid';'LFroMidO';'RFroMidO';'LFroInfO';'RFroInfO';'LFroInfT'; ...
                'RFroInfT';'LFroInfO';'RFroInfO';'LRolOpe';'RRolOpe';'LSupMotA';'RSupMotA'; ...
                'LOlfacto';'ROlfacto';'LFroSupM';'RFroSupM';'LFroMidO';'RFroMidO';'LRectus'; ...
                'RRectus';'LInsula';'RInsula';'LCinAnt';'RCinAnt';'LCinMid';'RCinMid';'LCinPos'; ...
                'RCinPos';'LHippoca';'RHippoca';'LParaHip';'RParaHip';'LAmygdal';'RAmygdal'; ...
                'LCalcari';'RCalcari';'LCuneus';'RCuneus';'LLingual';'RLingual';'LOccSup'; ...
                'ROccSup';'LOccMid';'ROccMid';'LOccInf';'ROccInf';'LFusifor';'RFusifor'; ...
                'LPostcen';'RPostcen';'LParSup';'RParSup';'LParInf';'RParInf';'LSupraMa'; ...
                'RSupraMa';'LAngular';'RAngular';'LPrecune';'RPrecune';'LParLob';'RParLob'; ...
                'LCaudate';'RCaudate';'LPutamen';'RPutamen';'LPallidu';'RPallidu';'LThalamu'; ...
                'RThalamu';'LHeschl';'RHeschl';'LTemSup';'RTemSup';'LTemPolS';'RTemPolS'; ...
                'LTemMid';'RTemMid';'LTemPolM';'RTemPolM';'LTemInf';'RTemInf';'LCerCru'; ...
                'RCerCru';'LCerCru';'RCerCru';'LCer3';'RCer3';'LCer45';'RCer45';'LCer6'; ...
                'RCer6';'LCer7b';'RCer7b';'LCer8';'RCer8';'LCer9';'RCer9';'LCer10';'RCer10'; ...
                'Ver12';'Ver3';'Ver45';'Ver6';'Ver7';'Ver8';'Ver9';'Ver10';};
            header.channels = char(Label);
            header.locs.labels = Label';
            header.ch_names = Labels;
        end
    end
end

if header.numchannels == 116 & strcmp(header.channels(1,:),'LPrecent')
    if ~exist('cfg','var') | ~isfield(cfg,'EXTRA') | ~isfield(cfg.EXTRA,'reduceAAL')
        if ~exist('cfg','var') | ~isfield(cfg,'MAIN') | ~isfield(cfg.MAIN,'auto') | cfg.MAIN.auto == 0
            cfg.EXTRA.reduceAAL = false;
            Prompt={'Reduce to 78 channels','reduceAAL'};
            Formats.type = 'check';
            cfg.EXTRA = inputsdlg(Prompt,'Reduce AAL',Formats,cfg.EXTRA);
        else
            cfg.EXTRA.reduceAAL = true;
        end
    end
    if cfg.EXTRA.reduceAAL == true
        [includechannels,includetext] = lab_get_AALselection(2);
        disp(['     atlas is AAL, reduce to standard regions (standard ' includetext ')'])
        [data,header] = lab_reduce_channels(data,header,includechannels);
        clearvars includechannels
    end
end

if isstruct(header) & isfield(header,'orig') & isfield(header.orig,'info') & isfield(header.orig.info,'dig')
    tmp = header.orig.info.dig;
    if ~isempty(tmp)
        header.locs.digits = zeros(size(tmp,2),3);
        for i = 1:size(tmp,2)
            header.locs.digits(i,:) = tmp(1,i).r' * 1000;
        end
    end
end

% lab_write_locs(fullfile(filepath,[filenameS '.els']),header.locs);

return






