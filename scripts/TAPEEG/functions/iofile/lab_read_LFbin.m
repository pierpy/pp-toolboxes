% LFbin = lab_read_LFbin(LFbin_file,numchannels,numsolutionpoints)
%
% Read Cartool Leadfields, number of channels and solutionpoints is needed
%
% written by F. Hatz 2012

function LFbin = lab_read_LFbin(LFbin_file,numchannels,numsolutionpoints)

if ~exist('LFbin_file','var') | ~ischar(LFbin_file)
    [LFbin_file,LFbin_filepath]=uigetfile('*.bin','Select LeadField file');
    [~,~,~,LFbin_fileS] = lab_filename(LFbin_file);
else
    [LFbin_file,LFbin_filepath,~,LFbin_fileS] = lab_filename(LFbin_file);
end
if strcmp(LFbin_fileS(end-2:end),'_LF')
    LFbin_fileS = LFbin_fileS(1:end-3);
end

if exist(fullfile(LFbin_filepath,[LFbin_fileS '.els']),'file')
    LOCS = lab_read_locs(fullfile(LFbin_filepath,[LFbin_fileS '.els']));
    if ~isempty(LOCS) & isfield(LOCS,'x')
        numchannels = length(LOCS.x);
    end
end
if ~exist('numchannels','var')
    numchannels = [];
end

if exist(fullfile(LFbin_filepath,[LFbin_fileS '.spi']),'file')
    SPI = lab_read_locs(fullfile(LFbin_filepath,[LFbin_fileS '.spi']));
    if ~isempty(SPI) & isfield(SPI,'x')
        numsolutionpoints = length(SPI.x);
    end
end
if ~exist('numsolutionpoints','var')
    numsolutionpoints = [];
end

if isempty(numchannels) | isempty(numsolutionpoints)
    prompt={'Number of channels','number of solutionpoints'};
    name='File characteristics';
    numlines(1:2,1) = ones(1,2);
    numlines(1:2,2) = ones(1,2)*25;
    defaultanswer={num2str(numchannels),num2str(numsolutionpoints)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    pause(0.1)
    numchannels = str2num(answer{1,1}); %#ok<ST2NM>
    numsolutionpoints = str2num(answer{2,1}); %#ok<ST2NM>
    clearvars 'name' 'answer' 'numlines' 'prompt' 'defaultanswer'
end

fid=fopen(fullfile(LFbin_filepath,LFbin_file));
for i = 1:numchannels
    tmp(:,i) = fread(fid,3 * numsolutionpoints,'float32');
end
fclose(fid);
tmp = reshape(tmp,3,numsolutionpoints,numchannels);
LFbin = permute(tmp,[3 1 2]);

return