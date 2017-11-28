function ROIS = lab_create_rois_shortlabel(ROIS)

if exist('ROIS','var') & ~isempty(ROIS) & ~isfield(ROIS,'label')
    for i = 1:ROIS.numrois
        tmp = ROIS.labels{i};
        if strfind(tmp,'_')
            tmp(strfind(tmp,'_')) = ' ';
            tmp = textscan(tmp,'%s');
            tmp = tmp{1,1};
            if strcmp(tmp{end,1}(1,1),'L') || strcmp(tmp{end,1}(1,1),'R')
                ROIS.label{i} = tmp{end,1}(1,1);
                tmp = tmp(1:end-1,1);
            else
                ROIS.label{i} = [];
            end
        else
            tmp = textscan(tmp,'%s');
            tmp = tmp{1,1};
            if size(tmp,1) > 1
                ROIS.label{i} = tmp{1,1}(1,1);
                tmp = tmp(2:end,1);
            else
                ROIS.label{i} = [];
            end
        end
        if size(tmp,1) > 1
            for j = 1:size(tmp,1)
                if length(tmp{j,1}) > 2
                    ROIS.label{i} = [ROIS.label{i} tmp{j,1}(1,1:3)];
                else
                    ROIS.label{i} = [ROIS.label{i} tmp{j,1}];
                end
            end
        else
            if length(tmp{1,1}) > 7
                ROIS.label{i} = [ROIS.label{i} tmp{1,1}(1,1:7)];
            else
                ROIS.label{i} = [ROIS.label{i} tmp{1,1}];
            end
        end
        if length(ROIS.label{i}) > 8
            ROIS.label{i} = ROIS.label{i}(1,1:8);
        end
    end
elseif ~exist('ROIS','var')
    ROIS = [];
end