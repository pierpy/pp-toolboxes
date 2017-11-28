% Plot results of BrainWave, permutation, brain regions (AAL-atlas) or connectivity matrix
% in source space, a file 'MNIbrain.mat' with atlas-template is needed
%
% lab_plot_IS(data,cfg)
%
% 'data' and 'cfg' are optional
%
% Written by F. Hatz Vumc Amsterdam 10/2012
% (distrubution only with permission of the author)

function lab_plot_IS(DATA,cfg)

if ~exist('cfg','var')
    cfg = [];
end
if ~exist('DATA','var')
    DATA = [];
else
    if isstruct(DATA) & isfield(DATA,'x')
        cfg.LOCS = DATA;
        DATA = [];
    elseif isnumeric(DATA)
        DATA2 = DATA;
        clearvars DATA
        if size(DATA2,1) > size(DATA2,2)
            DATA2 = DATA2';
        end
        if size(DATA2,1) == size(DATA2,2)
            for i = 1:size(DATA2,3)
                D1(i,:) = diag(DATA2(:,:,i)); %#ok<AGROW>
            end
            D2 = lab_extract_tril_wodiag(DATA2)';
            DATA = [];
            for i = 1:size(D2,1)
                DATA(end+1,1).data = D2(i,:); %#ok<AGROW>
                DATA(end,1).name = ['Data' num2str(length(DATA))];
                DATA(end,1).measure = ['Data' num2str(length(DATA))];
                DATA(end,1).subject = '';
                DATA(end,1).connections = true;
                DATA(end,1).nodes = false;
                DATA(end+1,1).data = D1(i,:); %#ok<AGROW>
                DATA(end,1).name = ['Data' num2str(length(DATA))];
                DATA(end,1).measure = ['Data' num2str(length(DATA))];
                DATA(end,1).subject = '';
                DATA(end,1).connections = false;
                DATA(end,1).nodes = nodes;
            end
        else
            for i = 1:size(DATA2,1)
                DATA(i,1).data = DATA2(1,:);
                DATA(i,1).name = ['Data' num2str(i)];
                DATA(i,1).measure = ['Data' num2str(i)];
                DATA(i,1).subject = '';
                DATA(end,1).connections = false;
                DATA(end,1).nodes = false;
            end
        end
        clearvars DATA2
    end
end

% set and plot data
cfg.DATA = DATA;
lab_set_plot_IS(cfg);

end

