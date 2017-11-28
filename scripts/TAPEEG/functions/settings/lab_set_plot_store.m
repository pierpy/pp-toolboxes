function settings = lab_set_plot_store(settings,matrix)

if ~exist('matrix','var')
    matrix = [];
end
if ~isfield(settings,'domatrix')
    settings.domatrix = 0;
end

if isfield(settings,'PLOT')
    PLOT = settings.PLOT;
end
if ~isfield(PLOT,'savepictures')
    PLOT(1).savepictures = false;
end
Prompt = cell(0,2);
Formats = [];
Prompt(end+1,:) = {'Save pictures' 'savepictures'};
Formats(end+1,1).type = 'check';
if settings.domatrix == 1
    if ~isfield(PLOT(1),'matrixmin1')
        if size(matrix,1) < 21 & length(unique(matrix(:))) > 2
            Prompt = cat(1,Prompt,{'Value 1 (default = all)','matrixmin1'; ...
                'Value 2 (default = 50%)','matrixmin2'; ...
                'Value 3 (default = 75%)','matrixmin3'});
            matrixtmp = lab_rm_diagonal(matrix);
            PLOT(1).matrixmin1 = 0;
            PLOT(1).matrixmin2 = median(matrixtmp(:));
            PLOT(1).matrixmin3 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin2));
            clearvars matrixtmp
        elseif length(unique(matrix(:))) > 2
            Prompt = cat(1,Prompt,{'Value 1 (default = 75%)','matrixmin1'; ...
                'Value 2 (default = 87.5%)','matrixmin2'; ...
                'Value 3 (default = 96.88%)','matrixmin3'});
            matrixtmp = lab_rm_diagonal(matrix);
            PLOT(1).matrixmin1 = median(matrixtmp(:));
            PLOT(1).matrixmin1 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin1));
            PLOT(1).matrixmin2 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin1));
            PLOT(1).matrixmin3 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin2));
            PLOT(1).matrixmin3 = median(matrixtmp(matrixtmp>PLOT(1).matrixmin3));
            clearvars matrixtmp
        else
            Prompt(end+1,:) = {'Value 1 (default = 0)','matrixmin1'};
            PLOT(1).matrixmin1 = 0;
            PLOT(1).matrixmin2 = [];
            PLOT(1).matrixmin3 = [];
        end
        if isnan(PLOT(1).matrixmin1)
            PLOT(1).matrixmin1 = 0;
        end
        if isnan(PLOT(1).matrixmin2)
            PLOT(1).matrixmin2 = [];
        end
        if isnan(PLOT(1).matrixmin3)
            PLOT(1).matrixmin3 = [];
        end
    else
        Prompt = cat(1,Prompt,{'Value 1','matrixmin1'; ...
            'Value 2','matrixmin2';'Value 3','matrixmin3'});
    end
    for i = 1:size(Prompt,1)-1
        Formats(end+1,1).type = 'edit'; %#ok<AGROW>
        Formats(end,1).format = 'float';
        Formats(end,1).limits = [0 inf];
    end
end
[PLOT(1),Cancelled] = inputsdlg(Prompt,'Save pictures',Formats,PLOT(1));
if Cancelled == 1
    PLOT(1).savepictures = false;
end
settings.PLOT = PLOT;

end