function eog = lab_get_eog(eog,header,cfg)

if exist('cfg','var') & isfield(cfg,'EXTRA') & isfield(cfg.EXTRA,'numdatachans')
    numchans = cfg.EXTRA.numdatachans;
elseif exist('header','var') & isfield(header,'numdatachannels')
    numchans = header.numdatachannels;
elseif exist('header','var') & isfield(header,'numchannels')
    numchans = header.numchannels;
else
    numchans = [];
end
if ~exist('eog','var') | isempty(eog)
    if numchans == 257
        eog = [37,241,46,244;18,238,10,234;31,252,32,253;31,226,25,225];
    elseif numchans == 214
        eog = [37,214,0,0;18,214,0,0];
    elseif numchans == 204
        eog = [37,204,0,0;18,204,0,0];
    else
        eog = [];
    end
end
if exist('header','var') & isfield(header,'channels')
    Chans = cellstr(header.channels)';
    Ceog = cell(size(eog));
    for i = 1:size(eog,1)
        for j = 1:size(eog,2)
            if ~isempty(eog(i,j)) & eog(i,j) > 0 & eog(i,j) <= length(Chans)
                Ceog{i,j} = Chans{eog(i,j)};
            else
                Ceog{i,j} = 'NaN';
            end
        end
    end
    Ceog = lab_table_dialog(Ceog,{'chan','ref','(chan)','(ref)'},'Channel list',1,{Chans,Chans,Chans,Chans});
    eog = zeros(size(Ceog));
    for i = 1:size(Ceog,1)
        for j = 1:size(Ceog,2)
            if ~isempty(Ceog{i,j}) & ~strcmp(Ceog{i,j},'NaN')
                eog(i,j) = find(strcmp(Chans,Ceog{i,j}));
            end
        end
    end
else
    eog = lab_table_dialog(eog,{'chan','ref','(chan)','(ref)'},'Channel list',1);
end
