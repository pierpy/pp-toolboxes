% Binarys matrix with different possible methods
% - fixed: fixed threshold
% - std: connecitivites > mean + threshold * std
% - percent: connecitivities stronger than threshold (in percent) 
%
% matrixout = lab_matrix2binary(matrix,threshold,mode)
%
% written by F. Hatz 2014

function matrixout = lab_matrix2binary(matrix,threshold,mode)

if ~exist('mode','var') | isempty(mode)
    mode = 'fixed';
end
if ~exist('threshold','var') | isempty(threshold)
    settings.mode = mode;
    settings.threshold = 0.5;
    Prompt = cell(0,2);
    Formats = [];
    
    Prompt(end+1,:) = {'Threshold','threshold'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'float';
    Formats(end,1).limits = [-inf inf];
    Formats(end,1).size = 40;
    
    Prompt(end+1,:) = {'','mode'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'fixed','std','percent'};
    Formats(end,1).callback = {@set_mode,'@ALL','@ALL'};
    
    [settings,Cancelled] = inputsdlg(Prompt,'Matrix2Binary',Formats,settings,2);
    if Cancelled == 1
        matrixout = matrix;
    end
    mode = settings.mode;
    threshold = settings.threshold;
end

switch mode
    case 'percent'
        for i = 1:size(matrix,3)
            tmp = lab_rm_diagonal(matrix(:,:,i));
            tmp = sort(tmp(:));
            tmp2 = round(length(tmp) * (threshold / 100));
            if tmp2 < 1
                tmp2 = 1;
            elseif tmp2 > length(tmp)
                tmp2 = length(tmp);
            end
            threshold2(i,1) = tmp(tmp2);
        end
        clearvars tmp tmp2 i
    case 'std'
        for i = 1:size(matrix,3)
            tmp = lab_rm_diagonal(matrix(:,:,i));
            threshold2(i,1) = mean(tmp(:)) + std(tmp(:)) * threshold;
        end
        clearvars tmp tmp2 i
    otherwise
        if size(threshold,2) == size(matrix,3)
            threshold2 = threshold';
        elseif size(threshold,1) ~= size(matrix,3)
            threshold2 = repmat(threshold(1),[size(matrix,3) 1]);
        else
            threshold2 = threshold;
        end
end

chans = size(matrix,1);
for i = 1:size(matrix,3)
    matrixtmp = matrix(:,:,i);
    matrixtmp2 = zeros(chans,chans);
    if threshold2(i,1) == 0
        matrixtmp2(matrixtmp>0) = 1;
    else
        matrixtmp2(matrixtmp>=threshold2(i,1)) = 1;
    end
    matrixout(:,:,i) = matrixtmp2;
end

end

function settings = set_mode(settings)
   if strcmp(settings.mode,'std')
       settings.threshold = 1.28;
   elseif strcmp(settings.mode,'percent')
       settings.threshold = 80;
   else
       settings.threshold = 0.5;
   end
end