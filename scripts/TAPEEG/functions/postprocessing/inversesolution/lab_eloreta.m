% Calculation of IS matrix (eLoreta) using fieldtrip code
%
% [ISmatrix,settings,SLOCS] = lab_eloreta(data,settings,header,SLOCS)
%
% Written by F. Hatz 2014

function [ISmatrix,settings,SLOCS] = lab_eloreta(data,settings,header,SLOCS)

if ~exist('SLOCS','var') | isempty(SLOCS)
    if exist('settings','var') & isfield(settings,'slocs')
        SLOCS = settings.slocs;
    else
        SLOCS = lab_load_spi;
    end
end

% Prepare Leadfield
if exist('settings','var') & isfield(settings,'LF') & isfield(settings.LF,'LFbin')
    Lp = settings.LF.LFbin;
else
    [settings.LF.LFbin_file,settings.LF.LFbin_filepath]=uigetfile('*.bin;*.mat;*.nii;*.hdr;*.fif;*.fiff','Select leadfield, mri or headmodel');
    if strcmp(settings.LF.LFbin_file(end-3:end),'.bin')
        if exist('SLOCS','var')
            numsolutionpoints = size(SLOCS.x,2);
        else
            prompt={'Number of solutionpoints'};
            name='File characteristics';
            numlines(1,1) = 1;
            numlines(1,2) = 25;
            answer=inputdlg(prompt,name,numlines);
            pause(0.1)
            numsolutionpoints = str2num(answer{1,1}); %#ok<ST2NM>
            clearvars 'name' 'answer' 'numlines' 'prompt'
        end
        Lp = lab_read_LFbin(fullfile(settings.LF.LFbin_filepath,settings.LF.LFbin_file),size(data,1),numsolutionpoints);
        settings.LF.LFbin = Lp;
    else
        if ~isempty(settings.LF.LFbin_file) & ~settings.LF.LFbin_file == 0
            settings.LF.mrifile = fullfile(settings.LF.LFbin_filepath,settings.LF.LFbin_file);
        end
        if ~isfield(settings.LF,'locs') & exist('header','var') & isfield(header,'locs')
            settings.LF.locs = header.locs;
        elseif ~isfield(settings.LF,'locs')
            settings.LF.locs = lab_read_locs;
        end
        [settings.LF] = lab_compute_leadfield(settings.LF,SLOCS,settings.LF.locs);
        Lp = settings.LF.LFbin;
    end
end

if ~isempty(Lp)
    % Normalize Leadfields
    disp('    Normalize leadfields')
    for i = 1:size(Lp,3)
        tmp = norm(Lp(:,:,i), 'fro');
        if tmp > 0
            Lp(:,:,i) = Lp(:,:,i) ./ tmp;
        end
    end
    
    % Collect variables for eLoreta
    if isfield(settings,'LF') & isfield(settings.LF,'LFbin')
        for i = 1:size(Lp,3)
            dipin.leadfield{1,i} = Lp(:,:,i);
        end
    end
    if isfield(settings,'LF') & isfield(settings.LF,'slocs')
        dipin.pos = settings.LF.slocs;
    else
        dipin.pos(:,1) = SLOCS.x';
        dipin.pos(:,2) = SLOCS.y';
        dipin.pos(:,3) = SLOCS.z';
    end
    dipin.inside = (1:size(dipin.pos,1))';
    dipin.outside = [];
    if ~isfield(settings,'LF') | ~isfield(settings.LF,'locs')
        if exist('header','var') & isfield(header,'locs')
            settings.LF.locs = header.locs;
        else
            settings.LF.locs = lab_read_locs;
        end
    end
    elec.chanpos(:,1) = settings.LF.locs.x';
    elec.chanpos(:,2) = settings.LF.locs.y';
    elec.chanpos(:,3) = settings.LF.locs.z';
    elec.elecpos = elec.chanpos;
    elec.label = settings.LF.locs.labels';
    elec.unit = 'mm';
    
    if isfield(settings,'LF') & isfield(settings.LF,'vol')
        vol = settings.LF.vol;
    else
        vol = [];
    end
    if ~exist('header','var') | ~isfield(header,'cov')
        header.cov = cov(data');
    end
    
    % Start eLoreta
    disp('    eLoreta (fieldtrip)')
    [dipout] = ft_eloreta(dipin,elec,vol,data,header.cov,'normalize','yes', ...
        'keepfilter','yes');

    % Collect results
    % if isfield(dipout,'leadfield')
    %     for i = 1:size(dipin.leadfield,2)
    %         settings.LF.LFbin(:,:,i) = dipin.leadfield{1,i};
    %     end
    % end
    
    if size(dipout.filter{1,1},1) == 1
        ISmatrix.matrix = zeros(size(dipout.filter,2),size(dipout.filter{1,1},2));
        for i = 1:size(dipout.filter,2)
            if ~isempty(dipout.filter{1,i})
                ISmatrix.matrix(i,:) = dipout.filter{1,i};
            end
        end
    elseif size(dipout.filter{1,1},1) == 3
        ISmatrix.x = zeros(size(dipout.filter,2),size(dipout.filter{1,1},2));
        ISmatrix.y = zeros(size(dipout.filter,2),size(dipout.filter{1,1},2));
        ISmatrix.z = zeros(size(dipout.filter,2),size(dipout.filter{1,1},2));
        for i = 1:size(dipout.filter,2)
            if ~isempty(dipout.filter{1,i})
                ISmatrix.x(i,:) = dipout.filter{1,i}(1,:);
                ISmatrix.y(i,:) = dipout.filter{1,i}(2,:);
                ISmatrix.z(i,:) = dipout.filter{1,i}(3,:);
            end
        end
    end
    ISmatrix.numsolutionpoints=size(dipout.filter,2);
    ISmatrix.numelectrodes=size(dipout.filter{1,1},2);
    ISmatrix.TSolutionPointName = cellstr(num2str((1:size(dipout.filter,2))'));
    ISmatrix.TElectrodeName = cellstr(num2str((1:size(dipout.filter{1,1},2))'));
    ISmatrix.regularization = 0;
    ISmatrix.type = 'eLoreta';
else
    ISmatrix = [];
end
