% Filtering EEG/MEG
%
% [data,header,settings,skip] = lab_filter(data,header,settings,novrb)
% [data,header,settings] = lab_filter(data,samplingrate,settings)
% [data,header,settings] = lab_filter(data,samplingrate)
% [data,header,settings] = lab_filter(data)
%
% data                 matrix (channels x samples)
% header               see 'lab_read_data'
% settings.filtermode  butter (butterworth)
%                      firls (high order least square)
%                      freq (frequency domain)
%                      wavelet (wavelet filter)
%                      cheby (Chebyshev Typ2)
% settings.lowpass     lowpass filter in Hz
% settings.highpass    highpass filter in Hz
% settings.notch       notch filter(s) in Hz
% settings.filtorder   Filter order (0 = automatic)
% novrb                if variable is defined, verbose-file is not written
%
% written by F. Hatz 2013


function [data,header,settings,skipprocessing] = lab_filter(data,header,settings,novrb)

skipprocessing = 0;

if ~exist('settings','var')
    settings = [];
end
if ~exist('header','var')
    tmp = inputdlg({'Samplingrate'},'Samplingrate',[1 25]);
    header.samplingrate = str2num(tmp{1}); %#ok<ST2NM>
    header.numchannels = size(data,1);
    header.numdatachannels = size(data,1);
    clearvars tmp
    novrb = 1; %#ok<NASGU>
elseif isnumeric(header)
    tmp.samplingrate = header;
    tmp.numchannels = size(data,1);
    tmp.numdatachannels = size(data,1);
    header = tmp;
    clearvars tmp
    novrb = 1; %#ok<NASGU>
end

if ~exist('novrb','var')
    if ~isfield(settings,'EEG_file') & isfield(header,'EEG_file')
        settings.EEG_file = header.EEG_file;
        settings.EEG_filepath = header.EEG_filepath;
    end
    if isfield(settings,'EEG_file')
        [~,~,~,settings.EEG_fileS] = lab_filename(settings.EEG_file);
    end
end
if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end

if ~isfield(settings,'filtermode')
    [settings,skipprocessing] = lab_set_filter(settings);
    if skipprocessing ~= 0
        return
    end
end
if strcmp(settings.filtermode,'none') | strcmp(settings.filtermode,'input')
    return
end
if ~isfield(settings,'nodisp')
    settings.nodisp = 0;
end
if ~isfield(settings,'filterauxchannels')
    settings.filterauxchannels = 0;
end

if ~isfield(settings,'notch') | settings.notch <= 0
    settings.notch = [];
end
if ~isfield(settings,'lowpass') | settings.lowpass <= 0
    settings.lowpass = [];
end
if ~isfield(settings,'highpass') | settings.highpass <= 0
    settings.highpass = [];
end
if settings.lowpass >= (header.samplingrate/2 - 2)
    settings.lowpass = (header.samplingrate/2 - 2);
    settings.lowpassA = (header.samplingrate/2 - 1);
elseif settings.lowpass > 0
    settings.lowpassA = settings.lowpass + 1;
end
if settings.highpass > 2
    settings.highpassA = settings.highpass - 1;
elseif settings.highpass > 0
    settings.highpassA = settings.highpass * 0.5;
end

if isfield(settings,'filterauxchannels') & settings.filterauxchannels == 1
    channels = header.numchannels;
elseif isfield(header,'numdatachannels')
    channels = header.numdatachannels;
end
if channels > size(data,1) | ~exist('channels','var')
    channels = size(data,1);
end

if (isempty(settings.lowpass) | isnan(settings.lowpass)) & ...
        (isempty(settings.highpass) | isnan(settings.highpass)) & isempty(settings.notch)
    return
end

for filternr = size(settings.notch,2):-1:0
    if settings.nodisp == 0 & ~isempty(settings.filtermode) & ~strcmp(settings.filtermode,'none')
        if filternr > 0
            disp(['   Filtering (Notch ' num2str(settings.notch(1,filternr)) 'Hz) - ' settings.filtermode])
        elseif isempty(settings.lowpass) | isnan(settings.lowpass)
            disp(['   Filtering (' num2str(settings.highpass) 'Hz highpass) - ' settings.filtermode])
        elseif isempty(settings.highpass) | isnan(settings.highpass)
            disp(['   Filtering (' num2str(settings.lowpass) 'Hz lowpass) - ' settings.filtermode])
        else
            disp(['   Filtering (' num2str(settings.highpass) 'Hz-' num2str(settings.lowpass) 'Hz) - ' settings.filtermode])
        end
    end
    
    % Butterworth filter
    if strcmp(settings.filtermode,'butter')
        if isfield(settings,'filtorder') & settings.filtorder > 0
            filtorder = round(settings.filtorder/2);
        else
            filtorder = 2;
        end
        
        warning off %#ok<WNOFF>
        testgood = 0;
        while testgood == 0
            if filternr > 0
                [B, A] = butter(filtorder,[(settings.notch(1,filternr)-2)/(header.samplingrate/2) (settings.notch(1,filternr)+2)/(header.samplingrate/2)],'stop');
            elseif isempty(settings.highpass) | isnan(settings.highpass)
                [B, A] = butter(filtorder,settings.lowpass/(header.samplingrate/2),'low');
            elseif isempty(settings.lowpass) | isnan(settings.lowpass)
                [B, A] = butter(filtorder,settings.highpass/(header.samplingrate/2),'high');
            else
                [B, A] = butter(filtorder,[settings.highpass/(header.samplingrate/2) settings.lowpass/(header.samplingrate/2)]);
            end
            settings.filtorder = filtorder*2;
            
            testgood = 1;
            tmp = filtfilt(B,A,data(1,:)-mean(data(1,:)));
            if find(isnan(tmp)) > 0
                testgood = 0;
                if filtorder > 20
                    filtorder = filtorder - 10;
                    disp('      filter error, reduce filter order by 10')
                else
                    filtorder = filtorder - 1;
                    disp('      filter error, reduce filter order by 1')
                end
            end
        end
        if settings.nodisp == 0
            disp(['      with order: ' num2str(filtorder*2)])
        end
        data(1,:) = tmp;
        warning on %#ok<WNON>
        for ch=2:channels
            data(ch,:)=filtfilt(B,A,data(ch,:)-mean(data(ch,:)));
        end
        header.filter.mode = 'butter';
        header.filter.order = filtorder*2;
        clearvars A B ch
        
    % Least square filter    
    elseif strcmp(settings.filtermode,'firls')
        Wstop1 = 1;     % First Stopband Weight
        Wpass  = 1;     % Passband Weight
        Wstop2 = 1;     % Second Stopband Weight
        if isfield(settings,'filtorder') & settings.filtorder > 0
            filtorder = round(settings.filtorder/2);
        else
            filtorder=round(2.4*header.samplingrate);
        end
        if size(data,2)/2 < filtorder*3;
            filtorder = floor(size(data,2)/6);
        end
        if rem(filtorder,2)~=0
            if filtorder > 2
                filtorder=filtorder-1;
            else
                filtorder=filtorder+1;
            end
        end
        
        warning off %#ok<WNOFF>
        testgood = 0;
        while testgood == 0
            % Calculate the coefficients using the FIRLS function.
            if filternr > 0
                hpfwts  = firls(filtorder, [0 (settings.notch(1,filternr)-3) (settings.notch(1,filternr)-2) (settings.notch(1,filternr)+2) (settings.notch(1,filternr)+3) header.samplingrate/2]/(header.samplingrate/2), [1 1 0 0 1 1], [Wstop1 Wpass Wstop2]);
            elseif isempty(settings.highpass) | isnan(settings.highpass)
                hpfwts  = firls(filtorder, [0 settings.lowpass settings.lowpassA header.samplingrate/2]/(header.samplingrate/2), [1 1 0 0], [Wpass Wstop2]);
            elseif isempty(settings.lowpass) | isnan(settings.lowpass)
                hpfwts  = firls(filtorder, [0 settings.highpassA settings.highpass header.samplingrate/2]/(header.samplingrate/2), [0 0 1 1], [Wstop1 Wpass]);
            else
                hpfwts  = firls(filtorder, [0 settings.highpassA settings.highpass settings.lowpass settings.lowpassA header.samplingrate/2]/(header.samplingrate/2), [0 0 1 1 0 0], [Wstop1 Wpass Wstop2]);
            end
            settings.filtorder = filtorder*2;
            
            testgood = 1;
            tmp = filtfilt(hpfwts,1,data(1,:)-mean(data(1,:)));
            if find(isnan(tmp)) > 0
                testgood = 0;
                if filtorder > 200
                    filtorder = filtorder - 100;
                    disp('      filter error, reduce filter order by 100')
                elseif filtorder > 20
                    filtorder = filtorder - 10;
                    disp('      filter error, reduce filter order by 10')
                else
                    filtorder = filtorder - 2;
                    disp('      filter error, reduce filter order by 2')
                end
            end
        end
        data(1,:) = tmp;
        if settings.nodisp == 0
            disp(['      with order: ' num2str(filtorder*2)])
        end
        warning on %#ok<WNON>
        for ch=2:channels
            data(ch,:)=filtfilt(hpfwts,1,data(ch,:)-mean(data(ch,:)));
        end
        clearvars hpfwts ch Wstop1 Wpass Wstop2;
        header.filter.mode = 'firls';
        header.filter.order = filtorder*2;
    
    % Frequency domain filter
    elseif strcmp(settings.filtermode,'freq')
        filter = zeros(1,size(data,2));
        freqbin = header.samplingrate / size(data,2);
        lowpassFD = round(settings.lowpass/freqbin);
        highpassFD = round(settings.highpass/freqbin);
        notchFD = round(settings.notch / freqbin);
        if filternr > 0
            filter = ones(1,size(data,2));
            filter(1,((notchFD(1,filternr)-2)+1:(notchFD(1,filternr)+2)+1)) = 0;
            filter(1,(end-(notchFD(1,filternr)+2)+1:end-(notchFD(1,filternr)-2)+1)) = 0;
        elseif isempty(settings.highpass) | isnan(settings.highpass)
            filter(1,(1:lowpassFD+1)) = 1;
            filter(1,(end-lowpassFD+1:end)) = 1;
        elseif isempty(settings.lowpass) | isnan(settings.lowpass)
            filter(1,(highpassFD+1:end-highpassFD+1)) = 1;
        elseif isfield(settings,'dofreqshift') & settings.dofreqshift == 1
            filter(1,(highpassFD:lowpassFD)) = 1; % Shifting of frequency window by 1 frequency bin
            filter(1,(end-lowpassFD:end-highpassFD)) = 1; % Shifting of frequency window by 1 frequency bin
        else
            filter(1,(highpassFD+1:lowpassFD+1)) = 1;
            filter(1,(end-lowpassFD+1:end-highpassFD+1)) = 1;
        end
        for ch=1:channels
            datafft=fft(data(ch,:)-mean(data(ch,:)));
            data(ch,:) = real(ifft(datafft .* filter));
        end
        clearvars filter datafft freqbin lowpassFD highpassFD notchFD
        header.filter.mode = 'freq';
    
    % wavelet filter
    elseif strcmp(settings.filtermode,'wavelet')
        fs = header.samplingrate;
        if ~isfield(settings,'wavsmoothing')
            settings.wavsmoothing = 4;
        end
        wavsmoothing = 2^nextpow2(settings.wavsmoothing);
        minfreq = 10 / fs * wavsmoothing; % minfreq at fs:500 & wavsmoothing:4 = 0.08
        maxfreq = fs / 2;
        scaledef.s0 = (centfrq('morl') * fs) / maxfreq;
        scaledef.ds = log2(centfrq('morl') / (scaledef.s0 * minfreq / fs))/(fs/wavsmoothing-1);
        scaledef.nb = round(fs/wavsmoothing);
        scaledef.type = 'pow';
        for ch=1:channels
            % sig.val = data(ch,:);
            % sig.period = 1/header.samplingrate;
            % datawv = cwtft(sig,'scales',scaledef);
            datawv = cwtft(data(ch,:),'scales',scaledef);
            if ~exist('filterwv','var')
                wvfreq = scal2frq(datawv.scales,'morl',1/header.samplingrate);
                if filternr > 0
                    tmp = abs(wvfreq - settings.notch(1,filternr));
                    filterwv = find(tmp == min(tmp));
                elseif isempty(settings.highpass) | isnan(settings.highpass)
                    filterwv = 1:find(wvfreq > settings.lowpass, 1, 'last' );
                elseif isempty(settings.lowpass) | isnan(settings.lowpass)
                    filterwv = find(wvfreq < settings.highpass, 1):size(wvfreq,2);
                else
                    filterwv = [1:find(wvfreq > settings.lowpass, 1, 'last') ...
                        find(wvfreq < settings.highpass, 1):size(wvfreq,2)];
                end
            end
            datawv.cfs(filterwv,:) = 0;
            %tmp = setdiff(1:size(datawv.scales,2),filterwv);
            %datawv.cfs = datawv.cfs(tmp,:);
            %datawv.scales = datawv.scales(1,tmp);
            data(ch,:) = icwtft(datawv);
            clearvars datawv
        end
        clearvars filterwv fs wavsmoothing
        header.filter.mode = 'wavelet';
        header.filter.smoothingfactor = settings.wavsmoothing;
        header.filter.pseudofreqs = wvfreq;
        clearvars wvfreq
        
    % Chebyshev Typ2 filter    
    elseif strcmp(settings.filtermode,'cheby')
        fs = header.samplingrate/2;
        if isfield(settings,'filtorder') & settings.filtorder > 0
            filtorder = round(settings.filtorder/2);
            dbstop = 40;
        else
            dbstop = 40;
            filtorder = 6;
            if filternr > 0
                while filtorder > 5
                    dbstop = dbstop - 5;
                    if dbstop > 0
                        filtorder = round(cheb2ord([(settings.notch(1,filternr)-2)/fs (settings.notch(1,filternr)+2)/fs], ...
                            [(settings.notch(1,filternr)-3)/fs (settings.notch(1,filternr)+3)/fs],1,dbstop)/2);
                    else
                        filtorder = 5;
                    end
                end
            elseif isempty(settings.highpass) | isnan(settings.highpass)
                while filtorder > 5
                    dbstop = dbstop - 5;
                    if dbstop > 0
                        filtorder = round(cheb2ord(settings.lowpass/fs,settings.lowpassA/fs,1,dbstop)/2);
                    else
                        filtorder = 5;
                    end
                end
            elseif isempty(settings.lowpass) | isnan(settings.lowpass)
                while filtorder > 5
                    dbstop = dbstop - 5;
                    if dbstop > 0
                        filtorder = round(cheb2ord(settings.highpass/fs,settings.highpassA/fs,1,dbstop)/2);
                    else
                        filtorder = 5;
                    end
                end
            else
                while filtorder > 5
                    dbstop = dbstop - 5;
                    if dbstop > 0
                        filtorder = round(cheb2ord([settings.highpass/fs settings.lowpass/fs], ...
                            [settings.highpassA/fs settings.lowpassA/fs],1,dbstop)/2);
                    else
                        filtorder = 5;
                    end
                end
            end
        end
        
        warning off %#ok<WNOFF>
        testgood = 0;
        while testgood == 0
            if filternr > 0
                [B, A] = cheby2(filtorder,80,[(settings.notch(1,filternr)-2)/fs (settings.notch(1,filternr)+2)/fs],'stop');
            elseif isempty(settings.highpass) | isnan(settings.highpass)
                [B, A] = cheby2(filtorder,dbstop,settings.lowpass/fs,'low');
            elseif isempty(settings.lowpass) | isnan(settings.lowpass)
                [B, A] = cheby2(filtorder,dbstop,settings.highpass/fs,'high');
            else
                [B, A] = cheby2(filtorder,dbstop,[settings.highpass/fs settings.lowpass/fs]);
            end
            settings.filtorder = filtorder*2;
            
            testgood = 1;
            tmp = filtfilt(B,A,data(1,:)-mean(data(1,:)));
            if find(isnan(tmp)) > 0
                testgood = 0;
                if filtorder > 20
                    filtorder = filtorder - 10;
                    disp('      filter error, reduce filter order by 10')
                else
                    filtorder = filtorder - 1;
                    disp('      filter error, reduce filter order by 1')
                end
            end
        end
        data(1,:) = tmp;
        if settings.nodisp == 0
            disp(['      with order: ' num2str(filtorder*2)])
        end
        warning on %#ok<WNON>
        for ch=2:channels
            data(ch,:)=filtfilt(B,A,data(ch,:)-mean(data(ch,:)));
        end
        data(isnan(data)) = 0;
        header.filter.order = filtorder*2;
        header.filter.mode = 'chevy';
        header.filter.dBstop = dbstop;
        clearvars A B ch
        
    % no valid filter mode found
    else
        if ~isempty(settings.filtermode) & ~strcmp(settings.filtermode,'none') & ~strcmp(settings.filtermode,'input')
            disp('   no valid filtermode selected')
            skipprocessing = 1;
        end
        return
    end
    
    % Write header info
    if skipprocessing == 0
        if filternr > 0
            if ~isfield(header,'notch')
                header.notch = [];
            end
            header.notch = [header.notch settings.notch(1,filternr)];
        else
            header.lowpass = settings.lowpass;
            header.highpass = settings.highpass;
        end
    end
end

% Write verbose file
if ~exist('novrb','var') & isfield(settings,'EEG_fileS')
    Verbose_file=[settings.EEG_fileS '_filt.vrb'];
    fid=fopen(fullfile(settings.EEG_filepath,Verbose_file),'w');
    fprintf(fid,'Filter EEG\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n');
    fprintf(fid,['EEG input file: ' settings.EEG_file]);
    fprintf(fid,'\n\n');
    fprintf(fid,['Filtermode: ' settings.filtermode]);
    if strcmp(settings.filtermode,'firls') | strcmp(settings.filtermode,'butter') | strcmp(settings.filtermode,'cheby')
        fprintf(fid,[' (order: ' num2str(filtorder*2) ')']);
    elseif strcmp(settings.filtermode,'wavelet')
        fprintf(fid,[' (smoothing factor: ' num2str(settings.wavsmoothing) ')']);
    end
    fprintf(fid,'\n\n');
    fprintf(fid,'Notch Filter: ');
    for i = 1:size(settings.notch,2)
        fprintf(fid,[num2str(settings.notch(1,i)) 'Hz ']);
    end
    fprintf(fid,'\n\n');
    fprintf(fid,['Hghpass Filter: ' num2str(settings.highpass) 'Hz']);
    fprintf(fid,'\n\n');
    fprintf(fid,['Lowpass Filter: ' num2str(settings.lowpass) 'Hz']);
    fprintf(fid,'\n\n');
    fprintf(fid,['Filter auxillary channels: '  num2str(settings.filterauxchannels)]);
    fprintf(fid,'\n');
    fclose(fid);
end

return