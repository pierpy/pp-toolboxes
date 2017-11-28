% Function to correct results of connectivity analysis in source space for
% false negativ connectivites due to adjacent sources
%
% Result = lab_calculate_connect_correctdistance(Result,header,cfg)
%
%   cfg.CONNECT.correctdistance   1 = enabled
%   header.locs                   structure with location of sources (lab_read_locs)
%
% Written by F. Hatz 2014 Neurology Basel

function Result = lab_calculate_connect_correctdistance(Result,header,cfg)

disp('    Correct matrices for distances')

if isempty(Result) | isempty(cfg) | ~isfield(cfg.CONNECT,'correctdistance') | cfg.CONNECT.correctdistance <= 0
    return
end

frequencies = fieldnames(Result);
tmp = [];
for i = 1:size(frequencies,1)
    if isstruct(Result.(frequencies{i,1})) & strcmp(frequencies{i,1}(1),'F')
        tmp = [tmp i];
    end
end
if isempty(tmp)
    frequencies = 0;
else
    frequencies = frequencies(tmp,1);
end

for Nfreq = 1:size(frequencies,1)
    if ~isnumeric(frequencies)
        result = Result.(frequencies{Nfreq,1});
    else
        result = Result;
    end
    
    if ~exist('distance','var')
        if isfield(header,'locs') & isfield(header.locs,'x') & ...
                size(header.locs.x,2) == header.numdatachannels
            distance = lab_distance(header);
            distance(distance >= cfg.CONNECT.correctdistance) = 0;
            distance(distance > 0) = 1;
        elseif isfield(result,'locs') & isfield(result.locs,'x') & ...
                size(result.locs.x,2) == header.numdatachannels
            distance = lab_distance(result);
            distance(distance >= cfg.CONNECT.correctdistance) = 0;
            distance(distance > 0) = 1;
        else
            return
        end
    end
    
    result = lab_correct_distance(result,distance);
    
    if ~isnumeric(frequencies)
        Result.(frequencies{Nfreq,1}) = result;
    else
        Result = result;
    end
end

end