function [hdr] = read_header_mff(filename)
    orig = [];
    
    % get header info from .bin files
    binfiles = dir(fullfile(filename, 'signal*.bin'));
    if isempty(binfiles)
      error('could not find any signal.bin in mff directory')
    end
    
    for iSig = 1:length(binfiles)
      signalname = binfiles(iSig).name;
      fullsignalname = fullfile(filename, signalname);
      orig.signal(iSig).blockhdr = read_mff_bin(fullsignalname);
    end
    
    % get hdr info from xml files
    ft_hastoolbox('XML4MAT', 1, 0);
    warning('off', 'MATLAB:REGEXP:deprecated') % due to some small code xml2struct
    xmlfiles = dir( fullfile(filename, '*.xml'));
    % disp('reading xml files to obtain header info...')
    for i = 1:numel(xmlfiles)
      if strcmpi(xmlfiles(i).name(1:2), '._') % Mac sometimes creates this useless files, don't use them
      elseif strcmpi(xmlfiles(i).name(1:6), 'Events') % don't read in events here, can take a lot of time, and we can do that in ft_read_event
      else
        fieldname = xmlfiles(i).name(1:end-4);
        filename_xml  = fullfile(filename, xmlfiles(i).name);
        orig.xml.(fieldname) = xml2struct(filename_xml);
      end
    end
    warning('on', 'MATLAB:REGEXP:deprecated')
    
    %make hdr according to FieldTrip rules
    hdr = [];
    Fs = zeros(length(orig.signal),1);
    nChans = zeros(length(orig.signal),1);
    nSamples = zeros(length(orig.signal),1);
    for iSig = 1:length(orig.signal)
      Fs(iSig)      = orig.signal(iSig).blockhdr(1).fsample(1);
      nChans(iSig)  = orig.signal(iSig).blockhdr(1).nsignals;
      % the number of samples per block can be different
      nSamples_Block = zeros(length(orig.signal(iSig).blockhdr),1);
      for iBlock  = 1:length(orig.signal(iSig).blockhdr)
        nSamples_Block(iBlock) = orig.signal(iSig).blockhdr(iBlock).nsamples(1);
      end
      nSamples(iSig) = sum(nSamples_Block);
    end
    
    if length(unique(Fs)) > 1 || length(unique(nSamples)) > 1
      error('Fs and nSamples should be the same in all signals')
    end
    
    hdr.Fs          = Fs(1);
    hdr.nChans      = sum(nChans);
    hdr.nSamplesPre = 0;
    hdr.nSamples    = nSamples(1);
    hdr.nTrials     = 1;
    
    %-get channel labels, otherwise create them
    if isfield(orig.xml, 'sensorLayout') % asuming that signal1 is hdEEG sensornet, and channels are in xml file sensorLayout
      for iSens = 1:numel(orig.xml.sensorLayout.sensors)
        if strcmp(orig.xml.sensorLayout.sensors(iSens).sensor.type, '0') %EEG chans
          % this should be consistent with ft_senslabel
          hdr.label{iSens} = ['E' num2str(orig.xml.sensorLayout.sensors(iSens).sensor.number)];
        elseif strcmp(orig.xml.sensorLayout.sensors(iSens).sensor.type, '1') %REF;
          hdr.label{iSens} = ['E' num2str(iSens)]; % to be consistent with other egi formats
        else
          %non interesting channels like place holders and COM
        end
      end
      %check if the amount of lables corresponds with nChannels in signal 1
      if length(hdr.label) == orig.signal(1).blockhdr(1).nsignals
        %good
      elseif length(hdr.label) > orig.signal(1).blockhdr(1).nsignals
        warning('found more lables in xml.sensorLayout than channels in signal 1, thus can not use info in sensorLayout, creating labels on the fly')
        for iSens = 1:orig.signal(1).blockhdr(1).nsignals
          % this should be consistent with ft_senslabel
          hdr.label{iSens} = ['E' num2str(iSens)];
        end
      else warning('found less lables in xml.sensorLayout than channels in signal 1, thus can not use info in sensorLayout, creating labels on the fly')
        for iSens = 1:orig.signal(1).blockhdr(1).nsignals
          % this should be consistent with ft_senslabel
          hdr.label{iSens} = ['E' num2str(iSens)];
        end
      end
      % get lables for other signals
      if length(orig.signal) == 2
        if isfield(orig.xml, 'pnsSet') % signal2 is PIB box, and lables are in xml file pnsSet
          nbEEGchan = length(hdr.label);
          for iSens = 1:numel(orig.xml.pnsSet.sensors)
            hdr.label{nbEEGchan+iSens} = num2str(orig.xml.pnsSet.sensors(iSens).sensor.name);
          end
          if length(hdr.label) == orig.signal(2).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
            %good
          elseif length(hdr.label) < orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
            % warning('found less lables in xml.pnsSet than channels in signal 2, labeling with s2_unknownN instead')
            for iSens = length(hdr.label)+1 : orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
              hdr.label{iSens} = ['s2_unknown', num2str(iSens)];
            end
          else warning('found more lables in xml.pnsSet than channels in signal 2, thus can not use info in pnsSet, and labeling with s2_eN instead')
            for iSens = orig.signal(1).blockhdr(1).nsignals+1 : orig.signal(1).blockhdr(1).nsignals + orig.signal(2).blockhdr(1).nsignals
              hdr.label{iSens} = ['s2_E', num2str(iSens)];
            end
          end
        else % signal2 is not PIBbox
          warning('creating channel labels for signal 2 on the fly')
          for iSens = 1:orig.signal(2).blockhdr(1).nsignals
            hdr.label{end+1} = ['s2_E', num2str(iSens)];
          end
        end
      elseif length(orig.signal) > 2
        % loop over signals and label channels accordingly
        warning('creating channel labels for signal 2 to signal N on the fly')
        for iSig = 2:length(orig.signal)
          for iSens = 1:orig.signal(iSig).blockhdr(1).nsignals
            if iSig == 1 && iSens == 1
              hdr.label{1} = ['s',num2str(iSig),'_E', num2str(iSens)];
            else
              hdr.label{end+1} = ['s',num2str(iSig),'_E', num2str(iSens)];
            end
          end
        end
      end
    else %no xml.sensorLayout present
      warning('no sensorLayout found in xml files, creating channel labels on the fly')
      for iSig = 1:length(orig.signal)
        for iSens = 1:orig.signal(iSig).blockhdr(1).nsignals
          if iSig == 1 && iSens == 1
            hdr.label{1} = ['s',num2str(iSig),'_E', num2str(iSens)];
          else
            hdr.label{end+1} = ['s',num2str(iSig),'_E', num2str(iSens)];
          end
        end
      end
    end
    
    % check if multiple epochs are present
    if isfield(orig.xml,'epochs') && length(orig.xml.epochs) > 1
      % add info to header about which sample correspond to which epochs, becasue this is quite hard for user to get...
      epochdef = zeros(length(orig.xml.epochs),3);
      for iEpoch = 1:length(orig.xml.epochs)
        if iEpoch == 1
          epochdef(iEpoch,1) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./1000./hdr.Fs)+1;
          epochdef(iEpoch,2) = round(str2double(orig.xml.epochs(iEpoch).epoch.endTime)./1000./hdr.Fs);
          epochdef(iEpoch,3) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./1000./hdr.Fs); %offset corresponds to timing
        else
          NbSampEpoch = round(str2double(orig.xml.epochs(iEpoch).epoch.endTime)./1000./hdr.Fs - str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./1000./hdr.Fs);
          epochdef(iEpoch,1) = epochdef(iEpoch-1,2) + 1;
          epochdef(iEpoch,2) = epochdef(iEpoch-1,2) + NbSampEpoch;
          epochdef(iEpoch,3) = round(str2double(orig.xml.epochs(iEpoch).epoch.beginTime)./1000./hdr.Fs); %offset corresponds to timing
        end
      end
      %warning('the data contains multiple epochs with possibly discontinuous boundaries. Added ''epochdef'' to hdr.orig defining begin and end sample of each epoch. See hdr.orig.xml.epochs for epoch details, use ft_read_header to obtain header or look in data.dhr.')
      % sanity check
      if epochdef(end,2) ~= hdr.nSamples
        error('number of samples in all epochs do not add up to total number of samples')
        %hdr.nSamples = epochdef(end,2);
      end
      orig.epochdef = epochdef;
    end
    hdr.orig      = orig;