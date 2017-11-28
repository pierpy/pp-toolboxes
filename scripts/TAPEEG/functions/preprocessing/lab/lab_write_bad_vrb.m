function lab_write_bad_vrb(Filename,settings,bad)

if isnumeric(Filename)
    fid = Filename;
    doclose = false;
else
    fid=fopen(Filename,'w');
    fprintf(fid,'Detection of bad channels\n');
    fprintf(fid,datestr(now,0));
    fprintf(fid,'\n\n');
    doclose = true;
end
    
if ~exist('bad','var')
    bad = [];
end
if isfield(bad,'filetxt')
    fprintf(fid,['Bad channels by file (' bad.filetxt '):\n']);
    fprintf(fid,num2str(bad.file));
else
    fprintf(fid,'Bad channels by file:\n');
    fprintf(fid,'disabled');
end
fprintf(fid,'\n\n');
if isfield(settings,'freqlim50') & settings.freqlim50 > 0
    fprintf(fid,['Bad Spectra 50Hz (Limit ' num2str(round(settings.freqlim50)) 'percent):\n']);
    if isfield(bad,'spectral50')
        fprintf(fid,num2str(bad.spectral50));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'freqlim60') & settings.freqlim60 > 0
    fprintf(fid,['Bad Spectra 60Hz (Limit ' num2str(round(settings.freqlim60)) 'percent):\n']);
    if isfield(bad,'spectral60')
        fprintf(fid,num2str(bad.spectral60));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'freqlimhigh') & settings.freqlimhigh > 0
    fprintf(fid,['Bad Spectra '  num2str(settings.spectshigh(1,1)) '-' num2str(settings.spectshigh(1,2)) ...
        'Hz (Limit ' num2str(round(settings.freqlimhigh)) 'percent):\n']);
    if isfield(bad,'spectralhigh')
        fprintf(fid,num2str(bad.spectralhigh));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'freqlimlow') & settings.freqlimlow > 0
    fprintf(fid,['Bad Spectra '  num2str(settings.spectslow(1,1)) '-' num2str(settings.spectslow(1,2)) ...
        'Hz (Limit ' num2str(round(settings.freqlimlow)) 'percent):\n']);
    if isfield(bad,'spectrallow')
        fprintf(fid,num2str(bad.spectrallow));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvaluebroken') & settings.zvaluebroken > 0
    fprintf(fid,['Bad broken (zValue: '  num2str(settings.zvaluebroken) '):\n']);
    if isfield(bad,'broken')
        fprintf(fid,num2str(bad.broken));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvaluevars') & settings.zvaluevars > 0
    fprintf(fid,['Bad variance (zValue: '  num2str(settings.zvaluevars) '):\n']);
    if isfield(bad,'variance')
        fprintf(fid,num2str(bad.variance));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvaluehurst') & settings.zvaluehurst > 0
    fprintf(fid,['Bad hurst (zValue: '  num2str(settings.zvaluehurst) '):\n']);
    if isfield(bad,'hurst')
        fprintf(fid,num2str(bad.hurst));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvaluemedian') & settings.zvaluemedian > 0
    fprintf(fid,['Bad median gradient (zValue: '  num2str(settings.zvaluemedian) '):\n']);
    if isfield(bad,'median')
        fprintf(fid,num2str(bad.median));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvalueamplitude') & settings.zvalueamplitude > 0
    fprintf(fid,['Bad amplitude (zValue: '  num2str(settings.zvalueamplitude) '):\n']);
    if isfield(bad,'amplitude')
        fprintf(fid,num2str(bad.amplitude));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvaluekurtosis') & settings.zvaluekurtosis > 0
    fprintf(fid,['Bad kurtosis (zValue: '  num2str(settings.zvaluekurtosis) '):\n']);
    if isfield(bad,'kurtosis')
        fprintf(fid,num2str(bad.kurtosis));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'zvaluecorr') & settings.zvaluecorr > 0
    fprintf(fid,['Bad topo correlation (zValue: '  num2str(settings.zvaluecorr) '):\n']);
    if isfield(bad,'corr')
        fprintf(fid,num2str(bad.corr));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'PEAK2MIN') & ~isempty(settings.PEAK2MIN)
    if strcmp(settings.PEAK2MIN.mode,'threshold')
        fprintf(fid,['Peak2min (' num2str(settings.PEAK2MIN.lowfreqpeak) '-' ...
            num2str(settings.PEAK2MIN.highfreqpeak) 'Hz - ' settings.PEAK2MIN.mode ...
            ' - factor:' num2str(settings.PEAK2MIN.factor) ' - min:' ...
            num2str(settings.PEAK2MIN.threshold) '):\n']);
        if isfield(bad,'peak2min')
            fprintf(fid,[num2str(sum(bad.peak2min==1)/sum(bad.peak2min~=1)) '%']);
        end
    else
        fprintf(fid,['Peak2min (' num2str(settings.PEAK2MIN.lowfreqpeak) '-' ...
            num2str(settings.PEAK2MIN.highfreqpeak) 'Hz - ' settings.PEAK2MIN.mode ...
            ' - factor:' num2str(settings.PEAK2MIN.factor) '):\n']);
        if isfield(bad,'peak2min')
            fprintf(fid,['Mean: ' num2str(mean(bad.peak2min))]);
        end
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'MicroCorr') & settings.MicroCorr == true
    fprintf(fid,'Microstates correlation:\n');
    if isfield(bad,'microcorr')
        fprintf(fid,num2str(mean(bad.microcorr)));
    end
end
if isfield(settings,'zvalueeye') & settings.zvalueeye > 0
    fprintf(fid,['Bad eye (zValue: '  num2str(settings.zvalueeye) '):\n']);
    if isfield(bad,'eye')
        fprintf(fid,num2str(bad.eye));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'ecgdetect') & ~isempty(settings.ecgdetect)
    if isfield(settings,'ecg_ch') & ~isempty(settings.ecg_ch)
        fprintf(fid,['Bad ecg (channel: '  num2str(settings.ecg_ch) '):\n']);
    elseif isfield(settings,'ecg_rate') & ~isempty(settings.ecg_rate)
        fprintf(fid,['Bad ecg (rate: '  num2str(settings.ecg_rate) '):\n']);
    end
    if isfield(bad,'QRS')
        fprintf(fid,num2str(bad.QRS));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'AVG') & ~isempty(settings.AVG) & isfield(settings,'AVGstd') & ~isempty(settings.AVGstd)
    fprintf(fid,['Bad std in average by marker (std = ' num2str(settings.AVGstd) '):\n']);
    if isfield(bad,'AVG')
        fprintf(fid,num2str(bad.AVG));
    end
    fprintf(fid,'\n\n');
end
if isfield(settings,'markerexclude') & ~isempty(settings.markerexclude)
    fprintf(fid,'Markers excluded:\n');
    fprintf(fid,'%s ',settings.markerexclude{:});
    fprintf(fid,'\n\n');
end
fprintf(fid,'Channels without signals:\n');
if isfield(bad,'nosignal')
    fprintf(fid,num2str(bad.nosignal));
end
if doclose == true
    fclose(fid);
end