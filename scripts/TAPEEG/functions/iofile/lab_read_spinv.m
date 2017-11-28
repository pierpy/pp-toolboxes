% Read sLoreta .spinv
%
% ISmatrix = lab_read_spinv(filename)
%
% written by F. Hatz 2012

function ISmatrix = lab_read_spinv(varargin)
if ~nargin
    [IS_file,IS_filepath]=uigetfile({'*.spinv'},'Select is file');
else
    IS_file=varargin{1};
    tmp=findstr(IS_file,filesep);
    if ~isempty(tmp)
        IS_filepath=IS_file(1:tmp(end));
        IS_file=IS_file(tmp(end)+1:end);
    else
        IS_filepath=[];
    end
    clearvars tmp;
end
fid=fopen(fullfile(IS_filepath,IS_file));
if fid > 0
    ISmatrix.numsolutionpoints=fread(fid,1,'float32');
    ISmatrix.numelectrodes=fread(fid,1,'float32');
    for i = 1:ISmatrix.numsolutionpoints
        ISmatrix.x(i,:) = fread(fid,ISmatrix.numelectrodes,'float32');
        ISmatrix.y(i,:) = fread(fid,ISmatrix.numelectrodes,'float32');
        ISmatrix.z(i,:) = fread(fid,ISmatrix.numelectrodes,'float32');
    end
    fclose(fid);
end
return