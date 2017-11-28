% Write Cartool .ris
%
% lab_write_ris(filename,RIS)
%
% written by F. Hatz 2012

function lab_write_ris(filename,RIS,header)

if ~exist('header','var')
    header = [];
end

[~,filepath,~,filenameS] = lab_filename(filename);
filename = [filenameS '.ris'];

if ~isstruct(RIS)
    tmp = RIS;
    clearvars RIS
    RIS.data = tmp;
    RIS.isscalar = 1;
    RIS.numsolutionpoints = size(tmp,1);
    RIS.numtimeframes = size(tmp,2);
    RIS.samplingrate = 0;
    clearvars tmp
elseif isfield(RIS,'x')
    RIS.isscalar = 0;
else
    RIS.isscalar = 1;
end
if isfield(RIS,'x') & ~isfield(RIS,'numsolutionpoints')
    RIS.numsolutionpoints = size(RIS.x,1);
    RIS.numtimeframes = size(RIS.x,2);
    RIS.samplingrate = 0;
end

RIS.version = 'RI01';
fid=fopen(fullfile(filepath,filename),'w');

fwrite(fid,RIS.version,'int8');
fwrite(fid,RIS.numsolutionpoints,'int32');
fwrite(fid,RIS.numtimeframes,'int32');
fwrite(fid,RIS.samplingrate,'float32');
fwrite(fid,RIS.isscalar,'int8');

if RIS.isscalar == 1
    fwrite(fid,RIS.data,'float32');
else
    tmp(:,:,1) = RIS.x;
    tmp(:,:,2) = RIS.y;
    tmp(:,:,3) = RIS.z;
    tmp = permute(tmp,[3 1 2]);
    fwrite(fid,tmp,'float32');
end
fclose(fid);

if isfield(RIS,'events')
    Marker_file=fullfile(filepath,[filename '.mrk']);
    lab_write_mrk(Marker_file,RIS);
end

if isfield(RIS,'locs')
    LOCS = RIS.locs;
elseif isfield(header,'locs')
    LOCS = header.locs;
end
if exist('LOCS','var')
    SPI_file = fullfile(filepath,filename);
    SPI_file = lab_write_spi(SPI_file,LOCS);
    if ~isempty(SPI_file) & isfield(RIS,'mrifile') & ~isempty(RIS.mrifile)
        [~,~,~,SPI_file] = lab_filename(SPI_file);
        fidout=fopen([filenameS '.lm'],'w');
        fprintf(fidout,[filename native2unicode([13 10])]);
        fprintf(fidout,[SPI_file '.spi' native2unicode([13 10])]);
        fprintf(fidout,[regexprep(RIS.mrifile,'\\','\\\') native2unicode([13 10])]);
        fclose(fidout);
    end
end

% Write Marker-file
if isfield(header,'events')
    Marker_file=fullfile(filepath,[filename '.mrk']);
    lab_write_mrk(Marker_file,header);
end

% Write EEGinfo-file (*.txt)
if isfield(header,'ref_chan') 
    if length(filename) <= 10 | ~strcmp(filename(end-10:end),'ICAtopo.sef')
        lab_write_eeginfo(fullfile(filepath,filename),header)
    end
end

% Write individual freqbands
if isfield(header,'IFREQ') & ~isempty(header.IFREQ)
    IFREQ = header.IFREQ; %#ok<NASGU>
    save(fullfile(filepath,[filenameS '.ifreq']),'IFREQ');
end

return