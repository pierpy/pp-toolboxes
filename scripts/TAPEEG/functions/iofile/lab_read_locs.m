% Read electrodes localications
%
% LOCS = lab_read_locs(filename)
%
% written by F. Hatz 2012

function LOCS = lab_read_locs(filename)

if ~exist('filename','var') | ~exist(filename,'file')
    [LOC_file,LOC_filepath]=uigetfile({'*.els;*.sfp;*.xyz;*.sxyz;*.spi'},'Select electrodes file');
    filename = fullfile(LOC_filepath,LOC_file);
end
[~,filenameS_path,format,filenameS] = lab_filename(filename);

if exist(filename,'file')
    if exist(['lab_read_' format]) == 2 %#ok<EXIST>
        disp(['    read electrodes locations (' filename ')'])
        % script for reading files is selected by name (all scripts in folder 'iofile')
        tmp = 1;
        eval(['tmp = nargout(@lab_read_' format ');']);
        if tmp > 1
            try
                eval(['tmp2 = nargout(@lab_read_' format ');']);
                if tmp2 == 3
                    eval(['[~,tmp] = lab_read_' format '(filename,[],1);']);
                else
                    eval(['[~,tmp] = lab_read_' format '(filename);']);
                end
                if isfield(tmp,'locs')
                    LOCS = tmp.locs;
                else
                    LOCS = [];
                    return
                end
                if isfield(tmp,'digits')
                    LOCS.digits = tmp.digits;
                end
            catch err
                disp(getReport(err))
                LOCS = [];
                return
            end
        else
            try
                eval(['[LOCS] = lab_read_' format '(filename);']);
            catch err
                disp(getReport(err))
                LOCS = [];
                return
            end
        end
        clearvars tmp
        
        if isempty(LOCS) | ~isfield(LOCS,'x')
            LOCS = [];
            disp('       format of loc file not supported')
            return
        end
        
        % check for mm scale
        if max(abs(LOCS.x)) < 20 & max(abs(LOCS.y)) < 20 & max(abs(LOCS.z)) < 20
            disp('       locs are in cm, correct to mm')
            LOCS.x = LOCS.x * 10;
            LOCS.y = LOCS.y * 10;
            LOCS.z = LOCS.z * 10;
        end
        
        % correct orientation to RAS
        if sum(LOCS.x)^2 > sum(LOCS.y)^2 & (strcmp(format,'els') | strcmp(format,'xyz') | strcmp(format,'sfp'))
            disp('       correct orientation to RAS')
            tmp = LOCS.x;
            LOCS.x = -LOCS.y;
            LOCS.y = tmp;
            clearvars tmp
            LOCS = lab_locs2sph(LOCS);
        end
        if ~isfield(LOCS,'aux')
            LOCS.aux = 0;
        end
        if exist(fullfile(filenameS_path,[filenameS '.grad']),'file')
            disp(['       read gradiometer / digitizer (' filenameS '.grad)'])
            load(fullfile(filenameS_path,[filenameS '.grad']),'-mat');
            if exist('grad','var')
                LOCS.grad = grad;
            end
            if exist('digits','var')
                LOCS.digits = digits;
            end
        end
    else
        disp('       format of loc file not supported')
        LOCS = [];
    end
else
    disp('       loc file not found')
    LOCS = [];
end

