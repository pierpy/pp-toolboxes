function events = lab_read_vmrk(filename)

events.POS = [];
events.DUR = [];
events.OFF = [];
events.TYP = [];

fid=fopen(filename,'rt');
if fid==-1,
    error('cannot open BrainVision marker file')
end
line = [];
while ischar(line) || isempty(line)
    line = fgetl(fid);
    if ~isempty(line) && ~(isnumeric(line) && line==-1)
        if strncmpi(line, 'Mk', 2)
            % this line contains a marker
            tok = tokenize(line, '=', 0);    % do not squeeze repetitions of the seperator
            if length(tok)~=2
                warning('skipping unexpected formatted line in BrainVision marker file');
            else
                % the line looks like "MkXXX=YYY", which is ok
                % the interesting part now is in the YYY, i.e. the second token
                tok = tokenize(tok{2}, ',', 0);    % do not squeeze repetitions of the seperator
                if isempty(tok{1})
                    tok{1}  = [];
                end
                if isempty(tok{2})
                    tok{2}  = [];
                end
                events.TYP = [events.TYP cellstr([tok{1} tok{2}])];
                events.POS = [events.POS int64(str2num(tok{3}))]; %#ok<ST2NM>
                events.DUR = [events.DUR int64(str2num(tok{4}))]; %#ok<ST2NM>
                events.OFF = [events.OFF int64(0)];
            end
        end
    end
end
fclose(fid);