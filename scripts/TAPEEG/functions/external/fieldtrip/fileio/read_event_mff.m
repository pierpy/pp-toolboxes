function [event] = read_event_mff(filename,hdr)
    if isempty(hdr)
      hdr = read_header_mff(filename);
    end
    
    % get event info from xml files
    ft_hastoolbox('XML4MAT', 1, 0);
    warning('off', 'MATLAB:REGEXP:deprecated') % due to some small code xml2struct
    xmlfiles = dir( fullfile(filename, '*.xml'));
    %disp('reading xml files to obtain event info... This might take a while if many events/triggers are present')
    for i = 1:numel(xmlfiles)
        if strcmpi(xmlfiles(i).name(1:6), 'Events')
            fieldname = xmlfiles(i).name(1:end-4);
            filename_xml  = fullfile(filename, xmlfiles(i).name);
            xml.(fieldname) = xml2struct(filename_xml);
        end
    end
    warning('on', 'MATLAB:REGEXP:deprecated')
    
    if exist('xml')
        % construct info needed for FieldTrip Event
        eventNames = fieldnames(xml);
        begTime = hdr.orig.xml.info.recordTime;
        begTime(11) = ' '; begTime(end-5:end) = [];
        begSDV = datenum(begTime);
        % find out if there are epochs in this dataset
        if isfield(hdr.orig.xml,'epochs') && length(hdr.orig.xml.epochs) > 1
            Msamp2offset = zeros(2,size(hdr.orig.epochdef,1),1+max(hdr.orig.epochdef(:,2)-hdr.orig.epochdef(:,1)));
            Msamp2offset(:) = NaN;
            for iEpoch = 1:size(hdr.orig.epochdef,1)
                nSampEpoch = hdr.orig.epochdef(iEpoch,2)-hdr.orig.epochdef(iEpoch,1)+1;
                Msamp2offset(1,iEpoch,1:nSampEpoch) = hdr.orig.epochdef(iEpoch,1):hdr.orig.epochdef(iEpoch,2); %sample number in samples
                Msamp2offset(2,iEpoch,1:nSampEpoch) = hdr.orig.epochdef(iEpoch,3):hdr.orig.epochdef(iEpoch,3)+nSampEpoch-1; %offset in samples
            end
        end
        
        % construct event according to FieldTrip rules
        eventCount = 0;
        for iXml = 1:length(eventNames)
            for iEvent = 1:length(xml.(eventNames{iXml}))
                if isfield(xml.(eventNames{iXml})(iEvent),'event')
                    eventTime  = xml.(eventNames{iXml})(iEvent).event.beginTime;
                    eventTime(11) = ' '; eventTime(end-5:end) = [];
                    if strcmp('-',eventTime(21))
                        % event out of range (before recording started): do nothing.
                    else
                        eventSDV = datenum(eventTime);
                        eventOffset = round((eventSDV - begSDV)*24*60*60*hdr.Fs); %in samples, relative to start of recording
                        if eventOffset < 0
                            % event out of range (before recording started): do nothing
                        else
                            eventCount = eventCount+1;
                            % calculate eventSample, relative to start of epoch
                            if isfield(hdr.orig.xml,'epochs') && length(hdr.orig.xml.epochs) > 1
                                for iEpoch = 1:size(hdr.orig.epochdef,1)
                                    [dum,dum2] = intersect(squeeze(Msamp2offset(2,iEpoch,:)), eventOffset);
                                    if ~isempty(dum2)
                                        EpochNum = iEpoch;
                                        SampIndex = dum2;
                                    end
                                end
                                eventSample = Msamp2offset(1,EpochNum,SampIndex);
                            else
                                eventSample = eventOffset+1;
                            end
                            
                            event(eventCount).type     = eventNames{iXml}(8:end);
                            event(eventCount).sample   = eventSample;
                            event(eventCount).offset   = 0;
                            event(eventCount).duration = str2double(xml.(eventNames{iXml})(iEvent).event.duration)./1000000000*hdr.Fs;
                            event(eventCount).value    = xml.(eventNames{iXml})(iEvent).event.code;
                        end  %if that takes care of non "-" events that are still out of range
                    end %if that takes care of "-" events, which are out of range
                end
            end %iEvent
        end
    else
        event = [];
    end