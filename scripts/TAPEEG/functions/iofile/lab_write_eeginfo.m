% Write TAPEEG eeginfo .txt
%
% lab_write_eeginfo(header,filename)
%
% written by F. Hatz 2012

function lab_write_eeginfo(filename,header)

if length(filename) > 7 & strcmp(filename(end-7:end),'_ICA.sef')
    isactivation = true;
else
    isactivation = false;
end
[~,filepath,~,filename] = lab_filename(filename);
filename = fullfile(filepath,[filename '.info']);

if ~isfield(header,'numdatachannels')
    header.numdatachannels = header.numchannels;
end
if ~isfield(header,'numauxchannels')
    header.numauxchannels = 0;
end

fid=fopen(filename,'w');
if isactivation == false
    fprintf(fid,'Bad channels');
else
    fprintf(fid,'Bad activations');
end
fprintf(fid,'\r\n');
if isfield(header,'badchans')
    fprintf(fid,num2str(header.badchans));
end
fprintf(fid,'\r\n');
if isactivation == false
    fprintf(fid,'Good channels');
else
    fprintf(fid,'Good activations');
end
fprintf(fid,'\r\n');
if isfield(header,'goodchans')
    fprintf(fid,num2str(header.goodchans));
else
    fprintf(fid,num2str(1:header.numdatachannels));
end
fprintf(fid,'\r\n');
if isfield(header,'interpolated')
    fprintf(fid,'Interpolated channels');
    fprintf(fid,'\r\n');
    fprintf(fid,num2str(header.interpolated));
    fprintf(fid,'\r\n');
end
fprintf(fid,'Reference channel');
fprintf(fid,'\r\n');
if isfield(header,'ref_chan') & isnumeric(header.ref_chan)
    fprintf(fid,num2str(header.ref_chan));
elseif isfield(header,'ref_chan')
    fprintf(fid,header.ref_chan);
end
fprintf(fid,'\r\n');
fprintf(fid,'Samplingrate');
fprintf(fid,'\r\n');
if isfield(header,'samplingrate')
    fprintf(fid,num2str(header.samplingrate));
end
fprintf(fid,'\r\n');
fprintf(fid,'Highpass filter');
fprintf(fid,'\r\n');
if isfield(header,'highpass')
    fprintf(fid,num2str(header.highpass));
end
fprintf(fid,'\r\n');
fprintf(fid,'Lowpass filter');
fprintf(fid,'\r\n');
if isfield(header,'lowpass')
    fprintf(fid,num2str(header.lowpass));
end
fprintf(fid,'\r\n');
fprintf(fid,'ECG channel');
fprintf(fid,'\r\n');
if isfield(header,'ecg_ch')
    fprintf(fid,num2str(header.ecg_ch));
end
fprintf(fid,'\r\n');
fprintf(fid,'Number EEG channels');
fprintf(fid,'\r\n');
fprintf(fid,num2str(header.numdatachannels));
fprintf(fid,'\r\n');
fprintf(fid,'Number AUX channels');
fprintf(fid,'\r\n');
fprintf(fid,num2str(header.numauxchannels));
fprintf(fid,'\r\n');
fprintf(fid,'Datatype');
fprintf(fid,'\r\n');
if isfield(header,'datatype')
    fprintf(fid,header.datatype);
else
    fprintf(fid,'');
end
if isfield(header,'quality')
    fprintf(fid,'\r\n');
    fprintf(fid,'Quality');
    fprintf(fid,'\r\n');
    fprintf(fid,num2str(header.quality));
end
if isfield(header,'begin')
    fprintf(fid,'\r\n');
    fprintf(fid,'Starttime in seconds');
    fprintf(fid,'\r\n');
    fprintf(fid,num2str(header.begin));
end
if isfield(header,'end')
    fprintf(fid,'\r\n');
    fprintf(fid,'Stoptime in seconds');
    fprintf(fid,'\r\n');
    fprintf(fid,num2str(header.end));
end
if isfield(header,'activationsexcluded')
    fprintf(fid,'\r\n');
    fprintf(fid,'Excluded activations');
    fprintf(fid,'\r\n');
    fprintf(fid,num2str(header.activationsexcluded));
end
fprintf(fid,'\r\n');
fprintf(fid,'Bad channels:');
fprintf(fid,'\r\n');
if isfield(header,'badchans') & ~isempty(header.badchans) & header.badchans > 0 & isfield(header,'locs')
    tmp = header.badchans(header.badchans<size(header.locs.labels,2));
    loctmp = header.locs.labels(1,tmp);
    for i = 1:size(tmp,2)
        fprintf(fid,[char(loctmp(1,i)) '\t']);
    end
else
    fprintf(fid,'');
end
fprintf(fid,'\r\n');
fclose(fid);