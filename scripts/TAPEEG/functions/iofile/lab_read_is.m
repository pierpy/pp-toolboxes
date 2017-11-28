% [ISmatrix,regularization] = lab_read_is(filename,regularization) - read Cartool *.is
%
% if variable 'regularization is set, only specified regularization is taken
%
% written by F. Hatz 2012

function [ISmatrix,regularization] = lab_read_is(filename,regularization)

if ~exist('regularization','var')
    flagall = 1;
    regularization = 0;
else
    flagall = 0;
end

if ~exist('filename','var') | ~exist(filename,'file')
    [filename,filepath]=uigetfile({'*.is;*.spinv'},'Select is file');
else
    tmp=strfind(filename,filesep);
    if ~isempty(tmp)
        filepath=filename(1:tmp(end));
        filename=filename(tmp(end)+1:end);
    else
        filepath=[];
    end
    clearvars tmp;
end

if strcmp(filename(end-2:end),'.is')
    fid=fopen(fullfile(filepath,filename));
    if fid > 0
        ISmatrix.version=strcat(fread(fid,4,'int8=>char')');
        if strcmp(ISmatrix.version,'IS03')
            ISmatrix.numelectrodes=fread(fid,1,'int32');
            ISmatrix.numsolutionpoints=fread(fid,1,'int32');
            ISmatrix.numregularizations=fread(fid,1,'int32');
            ISmatrix.isscalar=fread(fid,1,'int8');
            for i = 1:ISmatrix.numelectrodes
                ISmatrix.TElectrodeName(i,1)=cellstr(strcat(fread(fid,32,'int8=>char')'));
            end
            for i = 1:ISmatrix.numsolutionpoints
                ISmatrix.TSolutionPointName(i,1)=cellstr(strcat(fread(fid,16,'int8=>char')'));
            end
            for i = 1:ISmatrix.numregularizations
                ISmatrix.RegularizationValues(i,1)=fread(fid,1,'double');
            end
            for i = 1:ISmatrix.numregularizations
                ISmatrix.TRegularizationName(i,1)=cellstr(strcat(fread(fid,32,'int8=>char')'));
            end
            if flagall == 1
                regularization = ISmatrix.numregularizations-1;
            end
            if isempty(regularization) | regularization > (ISmatrix.numregularizations-1) | regularization == -1
                regularization = 0;
            end
            regnr = 1;
            for j = 1:(regularization+1)
                for i = 1:ISmatrix.numsolutionpoints
                    if ISmatrix.isscalar == 0
                        ISmatrix.x(i,:,regnr) = fread(fid,ISmatrix.numelectrodes,'float32');
                        ISmatrix.y(i,:,regnr) = fread(fid,ISmatrix.numelectrodes,'float32');
                        ISmatrix.z(i,:,regnr) = fread(fid,ISmatrix.numelectrodes,'float32');
                    else
                        ISmatrix.matrix(i,:,regnr) = fread(fid,ISmatrix.numelectrodes,'float32');
                    end
                end
                if flagall == 1
                    regnr = regnr + 1;
                end
            end
            if flagall == 0
                ISmatrix.regularization = regularization;
            end
            ISmatrix.type = 'Cartool';
        end
        fclose(fid);
    end
elseif strcmp(filename(end-4:end),'spinv')
    fid=fopen(fullfile(filepath,filename));
    if fid > 0
        ISmatrix.numsolutionpoints=fread(fid,1,'float32');
        ISmatrix.numelectrodes=fread(fid,1,'float32');
        for i = 1:ISmatrix.numsolutionpoints
            ISmatrix.x(i,:) = fread(fid,ISmatrix.numelectrodes,'float32');
            ISmatrix.y(i,:) = fread(fid,ISmatrix.numelectrodes,'float32');
            ISmatrix.z(i,:) = fread(fid,ISmatrix.numelectrodes,'float32');
        end
        fclose(fid);
        ISmatrix.TSolutionPointName = cellstr(num2str((1:ISmatrix.numsolutionpoints)'));
        ISmatrix.TElectrodeName = cellstr(num2str((1:ISmatrix.numelectrodes)'));
        ISmatrix.regularization = 0;
        ISmatrix.type = 'sLoreta';
        regularization = 0;
    end
else
    ISmatrix = [];
end
return