% Read marker files (.mrk,.edf,.xls)
%
% header = lab_read_markers(filename,header)
%
% filename = file to load
% header   = output of lab_read_data (if variable header is defined, marker
%            info stored in header will be combined with new marker info 
%            of loaded file)
%
% written by F. Hatz 2012

function [header,skipprocessing] = lab_read_markers(filename,header,novrb)

skipprocessing = 0;

if ~exist('novrb','var')
    novrb = false;
end
if ~exist('header','var')
    header.numtimeframes = 1000000000000;
end
[~,filepath,~,filenameS] = lab_filename(filename);
tmp = strfind(filenameS,'_');
if ~isempty(tmp) & tmp(1) ~= 1
    filenameS3 = filenameS(1:tmp(1)-1);
else
    filenameS3 = filenameS;
end
filenameS3 = fullfile(filepath,filenameS3);
filenameS2 = [fullfile(filepath,filenameS) filesep filenameS];
filenameS = fullfile(filepath,filenameS);
clearvars filepath

% Read .mrks files (Matlab-Container with markers)
markerfile = {[filename '.mrks'];[filenameS '.mrks'];[filenameS3 '.mrks']};
for i = 1:length(markerfile)
    if exist(markerfile{i},'file')
        if novrb == false
            disp(['    read marker from ' markerfile{i}])
        end
        headertmp = header;
        if isfield(headertmp,'events')
            headertmp = rmfield(headertmp,'events');
        end
        MRKS = load(markerfile{i},'-mat');
        header = lab_mix_markers(header,MRKS);
        header.events.mrkimport = 1;
        clearvars MRKS
        
    end
end

% Read .mrk files
markerfile = {[filename '.mrk'];[filename '.mrk1'];[filename '.mrk2']; ...
    [filenameS '.mrk'];[filenameS '.mrk1'];[filenameS '.mrk2']; ...
    [filenameS3 '.mrk'];[filenameS3 '.mrk1'];[filenameS3 '.mrk2']};
for i = 1:length(markerfile)
    if exist(markerfile{i},'file')
        if novrb == false
            disp(['    read marker from ' markerfile{i}])
        end
        headertmp = header;
        if isfield(headertmp,'events')
            headertmp = rmfield(headertmp,'events');
        end
        headertmp = lab_read_mrk(markerfile{i},headertmp);
        header = lab_mix_markers(header,headertmp);
        header.events.mrkimport = 1;
        clearvars headertmp
    end
end
clearvars markerfile

% Read .edf files
if ~strcmp(filename(end-3:end),'.edf')
    headertmp = [];
    if exist([filenameS '.edf'],'file')
        try
            [~,headertmp] = lab_read_edf([filenameS '.edf']);
        catch err
            if novrb == false
                disp(getReport(err))
            end
            skipprocessing = 1;
        end
    elseif exist([filenameS2 '.edf'],'file')
        try
            [~,headertmp] = lab_read_edf([filenameS2 '.edf']);
        catch err
            if novrb == false
                disp(getReport(err))
            end
            skipprocessing = 1;
        end   
    end
    if isfield(headertmp,'events')
        if novrb == false
            disp(['    replace markers by markers from ' filenameS '.edf'])
        end
        header.events = headertmp.events;
        header.events.edfimport = 1;
        clearvars headertmp
    end
end

% Read .xls files
if exist([filenameS '.xlsx'],'file')
    xlsfilename = [filenameS '.xlsx'];
elseif exist([filenameS '.xls'],'file')
    xlsfilename = [filenameS '.xls'];
elseif exist([filenameS2 '.xlsx'],'file')
    xlsfilename = [filenameS2 '.xlsx'];
elseif exist([filenameS2 '.xls'],'file')
    xlsfilename = [filenameS2 '.xls'];
end
if exist('xlsfilename','var')
    if ispc
        [~,~,xlstmp] = xlsread(xlsfilename);
    else
        [~,~,xlstmp] = xlsread(xlsfilename,1,'','basic');
    end
end
if exist('xlstmp','var')
    if strcmp(xlstmp{1},'Onset')
        if novrb == false
            disp(['    read marker from ' filenameS '.xls(x) file'])
        end
        for i = 2:size(xlstmp,1)
            headertmp.events.POS(1,i-1) = int64(xlstmp{i,1});
            if size(xlstmp,2) == 3
                tmp = int64(xlstmp{i,2});
                if tmp >= 1
                    headertmp.events.DUR(1,i-1) = tmp;
                else
                    headertmp.events.DUR(1,i-1) = int64(1);
                end
                headertmp.events.TYP(1,i-1) = xlstmp(i,3);
            else
                headertmp.events.DUR(1,i-1) = int64(1);
                headertmp.events.TYP(1,i-1) = xlstmp(i,2);
            end
            headertmp.events.OFF(1,i-1) = int64(0);
        end
        if max(headertmp.events.POS) <= header.numtimeframes
            header = lab_mix_markers(header,headertmp);
            header.events.xlsimport = 1;
        end
        clearvars headertmp
    end
    clearvars xlstmp
end

return
    