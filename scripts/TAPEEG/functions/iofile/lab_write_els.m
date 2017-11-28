% Write Cartool .els
%
% filename = lab_write_els(LOCS,filename,bad)
%
% if variable bad is set, els-file includes bad channels info
%
% written by F. Hatz 2012

function filename = lab_write_els(filename,LOCS,bad)

if ~exist('filename','var')
    [filename,filepath]=uiputfile('*.els','Select file and path');
    [~,~,~,filenameS] = lab_filename(filename);
else
    [~,filepath,~,filenameS] = lab_filename(filename);
end
if isempty(filepath)
    filepath = pwd;
end
filename = [filenameS '.els'];

if isfield(LOCS,'locs')
    % write header with electrodes structure
    if exist('bad','var') & strcmp(bad,'bad') & isfield(LOCS,'badchans') & ~isempty(LOCS.badchans) & LOCS.badchans ~= 0
        doBAD = 1;
    else
        doBAD = 0;
    end
    for i = 1:size(LOCS.locs.x,2)
        positions{i,1} = num2str(LOCS.locs.x(1,i),6); %#ok<AGROW>
        positions{i,2} = num2str(LOCS.locs.y(1,i),6); %#ok<AGROW>
        positions{i,3} = num2str(LOCS.locs.z(1,i),6); %#ok<AGROW>
        positions{i,4} = regexprep(LOCS.locs.labels{1,i},' ','_'); %#ok<AGROW>
        if  isfield(LOCS,'badchans') && ~isempty(LOCS.badchans) && max(LOCS.badchans == i) && doBAD == 1
            positions{i,5} = 'BAD'; %#ok<AGROW>
        else
            positions{i,5} = []; %#ok<AGROW>
        end
    end
    % look for auxillary channels
    if LOCS.numchannels > LOCS.numdatachannels
        clusters = 2;
    else
        clusters = 1;
    end
    fid = fopen(fullfile(filepath,filename),'wt');
    fprintf(fid,['ES01\n' int2str(LOCS.numchannels) '\n' int2str(clusters) '\n']);
    fprintf(fid,['Electrodes\n' int2str(size(positions,1)) '\n3\n']);
    for i = 1 : size(positions,1)
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',positions{i,:});
    end;    
    if LOCS.numchannels > LOCS.numdatachannels
        fprintf(fid,['Auxillary channels\n' int2str(LOCS.numchannels - LOCS.numdatachannels) '\n0\n']);
        for i = 1:(LOCS.numchannels - LOCS.numdatachannels)
            if isfield(LOCS.locs,'auxlabels') & length(LOCS.locs.auxlabels) >= i
                fprintf(fid,[int2str(0) '\t' int2str(0) '\t' int2str(i) '\t' LOCS.locs.auxlabels{i} '\n']);
            elseif isfield(LOCS,'channels')
                fprintf(fid,[int2str(0) '\t' int2str(0) '\t' int2str(i) '\t' regexprep(LOCS.channels(LOCS.numdatachannels+i,:),' ','') '\n']);
            else
                fprintf(fid,[int2str(0) '\t' int2str(0) '\t' int2str(i) '\tAUX' int2str(i) '\n']);
            end
        end
    end
    fclose(fid);
elseif isfield(LOCS,'x')
    % write TAPEEG/eeglab electrodes structure
    for i = 1:size(LOCS.x,2)
        positions{i,1} = num2str(LOCS.x(1,i),6); %#ok<AGROW>
        positions{i,2} = num2str(LOCS.y(1,i),6); %#ok<AGROW>
        positions{i,3} = num2str(LOCS.z(1,i),6); %#ok<AGROW>
        positions{i,4} = regexprep(LOCS.labels{1,i},' ','_'); %#ok<AGROW>
    end
    if isfield(LOCS,'aux') & LOCS.aux > 0
        clusters = 2;
        numchannels = size(positions,1) + LOCS.aux;
    else
        clusters = 1;
        numchannels = size(positions,1);
    end
    fid = fopen(fullfile(filepath,filename),'wt');
    fprintf(fid,['ES01\n' int2str(numchannels) '\n' int2str(clusters) '\n']);
    fprintf(fid,['Electrodes\n' int2str(size(positions,1)) '\n3\n']);
    for i = 1 : size(positions,1)
        fprintf(fid,'%s\t%s\t%s\t%s\t \n',positions{i,:});
    end;
    if clusters == 2
        fprintf(fid,['Auxillary channels\n' int2str(LOCS.aux) '\n0\n']);
        for i = 1:(LOCS.aux)
            if isfield(LOCS,'auxlabels') & length(LOCS.auxlabels) >= i
                fprintf(fid,[int2str(0) '\t' int2str(0) '\t' int2str(i) '\t' LOCS.auxlabels{i} '\n']);
            else
                fprintf(fid,[int2str(0) '\t' int2str(0) '\t' int2str(i) '\tAUX' int2str(i) '\n']);
            end
        end
    end
    fclose(fid);
elseif isfield(LOCS,'chanpos')
    % Write fieldtrip electrodes structure
    for i = 1:size(LOCS.chanpos,1)
        positions{i,1} = num2str(LOCS.chanpos(i,1),6); %#ok<AGROW>
        positions{i,2} = num2str(LOCS.chanpos(i,2),6); %#ok<AGROW>
        positions{i,3} = num2str(LOCS.chanpos(i,3),6); %#ok<AGROW>
        positions{i,4} = LOCS.label{i,1}; %#ok<AGROW>
    end
    fid = fopen(fullfile(filepath,filename),'wt');
    fprintf(fid,['ES01\n' int2str(size(positions,1)) '\n1\n']);
    fprintf(fid,['Electrodes\n' int2str(size(positions,1)) '\n3\n']);
    for i = 1 : size(positions,1)
        fprintf(fid,'%s\t%s\t%s\t%s\t \n',positions{i,:});
    end;
    fclose(fid);
elseif isnumeric(LOCS) & size(LOCS,2) == 3
    for i = 1:size(LOCS,1)
        positions{i,1} = num2str(LOCS(i,1),6); %#ok<AGROW>
        positions{i,2} = num2str(LOCS(i,2),6); %#ok<AGROW>
        positions{i,3} = num2str(LOCS(i,3),6); %#ok<AGROW>
        positions{i,4} = ['E' num2str(i,'%03d')]; %#ok<AGROW>
    end
    fid = fopen(fullfile(filepath,filename),'wt');
    fprintf(fid,['ES01\n' int2str(size(positions,1)) '\n1\n']);
    fprintf(fid,['Electrodes\n' int2str(size(positions,1)) '\n3\n']);
    for i = 1 : size(positions,1)
        fprintf(fid,'%s\t%s\t%s\t%s\t \n',positions{i,:});
    end;
    fclose(fid);
else
    filename = [];
end

return