% Function to read simbio mesh files 
% (original function crashes on windows 7)
%
% [Pnt,Hex,Labels] = lab_read_vista_mesh(Filename)
%
% written by F. Hatz 2014

function [Pnt,Hex,Labels] = lab_read_vista_mesh(Filename)

fid=fopen(Filename,'r');
inputline = fgetl(fid);
if ~strcmp(inputline(1:8),'V-data 2')
    Header = [];
    disp('Abort, unkown format')
end

Nr = 0;
while ~strcmp(inputline,'}')
    inputline = fgetl(fid);
    if strcmp(inputline,'	graph: graph {')
        Nr = Nr+1;
        Var = '';
    elseif strcmp(inputline(end),'{')
        Var = process_line(inputline);
    elseif strcmp(inputline(end),'}')
        Var = '';
    else
        [name,value] = process_line(inputline);
        if isempty(Var)
            Header(1,Nr).(name) = value;
        else
            Header(1,Nr).(Var).(name) = value;
        end
    end
end

StartPos = ftell(fid) + 2;

for j = 1:size(Header,2)
    Vars = fieldnames(Header(1,j));
    Vars = cat(1,cellstr('dat'),Vars);
    for i = 1:length(Vars);
        if strcmp(Vars{i},'dat')
            HeaderTmp = Header(1,j);
        else
            HeaderTmp = Header(1,j).(Vars{i});
        end
        if isstruct(HeaderTmp);
            skip = 0;
            if strcmp(HeaderTmp.repn,'ubyte')
                format = '*uint8';
                correct = 0;
                factor = 1;
            elseif strcmp(HeaderTmp.repn,'sbyte')
                format = '*uint8';
                correct = -128;
                factor = 1;
            elseif strcmp(HeaderTmp.repn,'float')
                format = '*float32';
                correct = 0;
                factor = 4;
            elseif strcmp(HeaderTmp.repn,'double')
                format = '*double';
                correct = 0;
                factor = 8;
            elseif strcmp(HeaderTmp.repn,'short')
                format = '*int16';
                correct = 0;
                factor = 2;
            elseif strcmp(HeaderTmp.repn,'long')
                format = '*int32';
                correct = 0;
                factor = 4;
            elseif strcmp(HeaderTmp.repn,'bit')
                format = '*bit1';
                correct = 0;
                factor = 1/8;
            else
                skip = 1;
            end
            if skip == 0
                fseek(fid,HeaderTmp.data+StartPos,'bof');
                dat = fread(fid,HeaderTmp.length/factor,format,0,'b');
                if isfield(HeaderTmp,'ncomponents')
                    ncolumns = HeaderTmp.ncomponents;
                elseif isfield(HeaderTmp,'nfields')
                    ncolumns = HeaderTmp.nfields + 2;
                else
                    ncolumns = 1;
                end
                if isfield(HeaderTmp,'ncolumns')
                    nrows = HeaderTmp.ncolumns;
                elseif isfield(HeaderTmp,'size')
                    nrows = HeaderTmp.size;
                else
                    nrows = 1;
                end
                Result(1,j).(Vars{i}) = double(reshape(dat(1:(ncolumns*nrows)),[ncolumns nrows ])' + correct);
            end
        end
    end
end

if size(Result,2) == 2 & isfield(Result,'dat') & ...
        size(Result(1,1).dat,2) == 7 & size(Result(1,2).dat,2) == 11
    Pnt = Result(1,1).dat(:,4:6);
    Hex = Result(1,2).dat(:,4:11);
    Labels = Result(1,2).matprops;
end

end


function [name,value] = process_line(tline)

tline = strtrim(tline);
tmp = strfind(tline,':');
if ~isempty(tmp) & tmp > 1 & tmp < length(tline)
    name = tline(1:tmp-1);
    name = regexprep(name,{'(',')','/','-'},'');
    value = tline(tmp+2:end);
    warning off %#ok<WNOFF>
    tmp = str2num(value);
    warning on %#ok<WNON>
    if ~isempty(tmp)
        value = tmp;
    end
else
    name = '';
    value = [];
end

end

