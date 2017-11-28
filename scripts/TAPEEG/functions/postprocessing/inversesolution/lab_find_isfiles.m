% Find files for claculation of inverse solution (template-files in
% standard / project folder)
%
% written by F. Hatz 2012

function settings = lab_find_isfiles(cfg,header,settings,doforward)

global Main_Path

if ~exist('doforward','var')
    doforward = false;
end
if ~exist('settings','var')
    settings = [];
end
if ~isfield(settings,'IS_file')
    settings.IS_file = '';
end
if ~isfield(settings,'SPI_file')
    settings.SPI_file = '';
end
if ~isfield(settings,'LOC_file')
    settings.LOC_file = '';
end
if ~isfield(settings,'MRI_file')
    settings.MRI_file = '';
end
if ~isfield(settings,'ROIS_file')
    settings.ROIS_file = '';
end

% find IS_file if empty
if isempty(settings.IS_file)
    if doforward == false
        strings = {'.hdr','.nii','iso.fif','iso.fiff','mri.fif','mri.fiff','.mat','LF.bin','.is','.spinv'};
    else
        strings = {'.hdr','.nii','iso.fif','iso.fiff','mri.fif','mri.fiff','.mat','LF.bin'};
    end
    if exist(fullfile(Main_Path,'IS'),'dir')
        searchpath = fullfile(Main_Path,'IS');
        settings.IS_file = searchfile(searchpath,strings);
    end
    if exist('header','var') & isfield(header,'numdatachannels') & exist([fullfile(Main_Path,'IS') num2str(header.numdatachannels)],'dir')
        searchpath = [fullfile(Main_Path,'IS') num2str(header.numdatachannels)];
        tmp = searchfile(searchpath,strings);
        if ~isempty(tmp)
            settings.IS_file = tmp;
        end
        clearvars tmp
    end
    searchpath = pwd;
    tmp = searchfile(searchpath,strings);
    if ~isempty(tmp)
        settings.IS_file = tmp;
    end
    clearvars tmp
    if exist('cfg','var') & isfield(cfg,'settings_path')
        searchpath = cfg.settings_path;
        tmp = searchfile(searchpath,strings);
        if ~isempty(tmp)
            settings.IS_file = tmp;
        end
        clearvars tmp
    end
    if exist('header','var') & isfield(header,'EEG_filepath')
        searchpath = header.EEG_filepath;
        tmp = searchfile(searchpath,strings);
        if ~isempty(tmp)
            settings.IS_file = tmp;
        end
        clearvars tmp
    end
    if exist('cfg','var') & isfield(cfg,'Output_filepath')
        searchpath = cfg.Output_filepath;
        tmp = searchfile(searchpath,strings);
        if ~isempty(tmp)
            settings.IS_file = tmp;
        end
        clearvars tmp
    end
end

if isempty(settings.IS_file)
    return
else
    [~,searchpath,ISformat] = lab_filename(settings.IS_file);
end

if isempty(settings.SPI_file)
    switch ISformat
        case {'is','spinv','mat','bin'}
            strings = {'.spi'};
        otherwise
            strings = {'.spi','.hdr','.nii','iso.fif','iso.fiff','mri.fif','mri.fiff'};
    end
    settings.SPI_file = searchfile(searchpath,strings);
end
[~,~,SPIformat] = lab_filename(settings.SPI_file);

if isempty(settings.LOC_file)
    strings = {'.els','.xyz'};
    settings.LOC_file = searchfile(searchpath,strings);
end

if isempty(settings.MRI_file)
    strings = {'.hdr','.nii','iso.fif','iso.fiff','mri.fif','mri.fiff'};
    settings.MRI_file = searchfile(searchpath,strings);
end

if isempty(settings.ROIS_file)
    switch SPIformat
        case {'hdr','nii','fif','fiff'}
            strings = {'.nii','.rois'};
        otherwise
            strings = {'.rois'};
    end
    settings.ROIS_file = searchfile(searchpath,strings);
end

return

function filename = searchfile(searchpath,strings)
   filename = '';
   for i = 1:size(strings,2)
       tmp = dir(fullfile(searchpath,['*' strings{1,i}]));
       if ~isempty(tmp)
           filename = fullfile(searchpath,tmp(1,1).name);
       end
   end
return