% Read Cartool .rois
%
% [ROIS] = lab_read_rois(Filename)
%
% filename = .rois file
%
% written by F. Hatz 2012

function ROIS = lab_read_rois(Filename)

if ~exist('Filename','var') | isempty(Filename)
    [ROIS_file,ROIS_filepath]=uigetfile('*.rois','Select ROIS-file');
    if isempty(ROIS_file) | ROIS_file == 0
        ROIS = [];
        return
    else
        Filename = fullfile(ROIS_filepath,ROIS_file);
    end
    clearvars ROIS_file ROIS_filepath
end

if ~exist(Filename,'file')
    ROIS = [];
    return
end

fid=fopen(Filename);
if fid > 0
    ROIS.version=textscan(fgetl(fid),'%s');
    ROIS.version = ROIS.version{1,1}{1,1};
    if strcmp(ROIS.version,'RO01')
        ROIS.numsolutionpts = str2num(fgetl(fid));
        ROIS.numrois = str2num(fgetl(fid));
        ROIS.solutionptsAll = [];
        for i = 1:ROIS.numrois
            ROIS.labels{i} = fgetl(fid);
            ROIS.solutionpts{i} = str2num(fgetl(fid));
            ROIS.solutionptsAll = [ROIS.solutionptsAll ROIS.solutionpts{i}];
            ROIS.solutionptsAll = sort(ROIS.solutionptsAll);
        end
    end
    fclose(fid);
end

% Create short label
ROIS = lab_create_rois_shortlabel(ROIS);