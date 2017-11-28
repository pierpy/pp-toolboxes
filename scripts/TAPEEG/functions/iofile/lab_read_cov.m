% Read covariance matrix
%
% cov=lab_read_cov(cov_file,cfg)
% 
% by F. Hatz, Neurology Basel

function cov = lab_read_cov(cov_file)

if exist('cov_file','var')
    [cov_file,cov_filepath] = lab_filename(cov_file);
    if isempty(cov_filepath)
        cov_filepath = pwd;
    end
else
    [cov_file,cov_filepath]=uigetfile('*.cov','Select covariance file');
    if cov_file == 0
        cov = [];
        return
    end
end

try
    cov = importdata(fullfile(cov_filepath,cov_file));
catch %#ok<CTCH>
    cov  = [];
end
   
if size(cov,1) ~= size(cov,2)
    cov = [];
end

return


