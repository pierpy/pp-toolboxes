% [data,header,cfg]=lab_read_asci(filename,cfg)
% written by F. Hatz 2012

function [data,header,cfg]=lab_read_asci(filename,cfg)

if ~exist('cfg','var')
    cfg = [];
end

if exist('filename','var')
    ASCI_file=filename;
    tmp=findstr(ASCI_file,filesep);
    if ~isempty(tmp)
        ASCI_filepath=ASCI_file(1:tmp(end));
        ASCI_file=ASCI_file(tmp(end)+1:end);
    else
        ASCI_filepath=pwd;
    end
else
    [ASCI_file,ASCI_filepath]=uigetfile('*.txt','Select ASCI_data.txt file');
end

dat1=importdata(fullfile(ASCI_filepath,ASCI_file));
dat2=importdata(fullfile(ASCI_filepath,[ASCI_file(1:end-8) 'signals.txt']));
data=dat1.data(2:end,:)';
header.channels=dat2.textdata(2:end,2);
header.samplingrate=dat2.data(1,1);
clearvars dat1 dat2

% Read annotations.txt
fid=fopen(fullfile(ASCI_filepath,[ASCI_file(1:end-8) 'annotations.txt']));
if fid > -1
  tline=fgetl(fid);
  tline=fgetl(fid);
  i = 1;
  while ~isnumeric(tline)
    if regexp(tline,'\t','once')
      tmp=regexp(tline,'\t');
    else
      tmp=regexp(tline,',');
    end
    header.events.POS(i,1)=int64(round(str2num(tline(1:tmp(1)-1))*samplingrate)+1);
    if tmp(2) > (tmp(1)+1)
      header.events.DUR(i,1)=int64(round(str2num(tline(tmp(1)+1:tmp(2)-1))*samplingrate)+1);
    else  
      header.events.DUR(i,1)=int64(1);
    end
    header.events.TYP(i,1)=cellstr(tline(tmp(2)+1:end));
    tline=fgetl(fid);
    i = i + 1;
  end
end
fclose(fid);
clearvars tline i tmp
