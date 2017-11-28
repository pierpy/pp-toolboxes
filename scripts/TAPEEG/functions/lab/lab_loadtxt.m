% loadtxt() - load ascii text file into numeric or cell arrays
%
% array = loadtxt(filename)
%
%  filename - name of the input file
%
%
% Author: Arnaud Delorme, CNL / Salk Institute, 29 March 2002
% Modified for TAPEEG: F. Hatz 2012

function array = lab_loadtxt(filename)

skipprocessing = 0;

% open the file
% -------------
if exist(filename,'file') ~= 2
    skipprocessing = 1;
    disp( ['file ' filename ' not found'] );
end

if skipprocessing == 0
    fid=fopen(filename,'r','ieee-le');
    if fid<0
        skipprocessing = 1;
        disp( ['file ' filename ' found but error while opening file'] );
    end
end
if skipprocessing == 0
    inputline = fgetl(fid);
    linenb = 1;
    while isempty(inputline) | inputline~=-1
        colnb = 1;
        if ~isempty(inputline)
            tabFirstpos = 1;
            while ~isempty(deblank(inputline))
                if tabFirstpos && length(inputline) > 1 && all(inputline(1) ~= [9 32]), tabFirstpos = 0; end;
                [tmp,inputline,tabFirstpos] = mystrtok(inputline, [9 32], tabFirstpos);
                tmp2 = str2double(tmp);
                if isnan( tmp2 )
                    array{linenb, colnb} = tmp;
                else
                    array{linenb, colnb} = tmp2;
                end
                colnb = colnb+1;
            end
            linenb = linenb +1;
        end
        inputline = fgetl(fid);
        if linenb > Inf
            inputline = -1;
        end
    end
    fclose(fid);
end

% problem strtok do not consider tabulation
% -----------------------------------------
function [str, strout, tabFirstpos] = mystrtok(strin, delim, tabFirstpos)
% remove extra spaces at the beginning
while any(strin(1) == delim) && strin(1) ~= 9 && strin(1) ~= ','
    strin = strin(2:end);
end;
% for tab and coma, consider empty cells
if length(strin) > 1 && any(strin(1) == delim)
    if tabFirstpos || any(strin(2) == delim)
        str = '';
        strout = strin(2:end);
        if strin(2) ~= 9 && strin(2) ~= ','
            tabFirstpos = 0;
            strout = strtrim(strout);
        end
    else
        [str, strout] = strtok(strin, delim);
    end
else
    [str, strout] = strtok(strin, delim);
end
