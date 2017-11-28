function [settings,skipprocessing] = lab_set_calculate_kmeans(settings,dosearch)

skipprocessing = 0;

if ~exist('dosearch','var')
    dosearch = 2;
end

if ~exist('settings','var') | isempty(settings)
    settings.doPCA = false;
    settings.maxclusters = 40;
    settings.WriteMatrices = false;
    settings.results = '';
end

Formats = {};
Prompt = cell(0,2);

if dosearch == 2
    settings.searchfolder = pwd;
    settings.searchmode = 'Matrix';
    settings.searchstring = '';
    settings.excludestring = '';
    Prompt(end+1,:) = {'Search Folder','searchfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'dir';
    
    Prompt(end+1,:) = {'Search ','searchmode'};
    Formats(end+1,1).type = 'list';
    Formats(end,1).style = 'popupmenu';
    Formats(end,1).format = 'input';
    Formats(end,1).items = {'strings','Matrix'};
    
    Prompt(end+1,:) = {'','searchstring'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).size = 300;
    
    Prompt{end+1,1} = '';
    Formats(end+1,1).type = 'text';
    
    Prompt{end+1,1} = 'Exclude strings';
    Formats(end+1,1).type = 'text';
    Prompt(end+1,:) = {'','excludestring'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'vector';
    Formats(end,1).size = 300;
    
    Prompt{end+1,1} = '';
    Formats(end+1,1).type = 'text';
elseif dosearch == 1
    Prompt(end+1,:) = {'Output-Folder','searchfolder'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'dir';
    
    Prompt{end+1,1} = '';
    Formats(end+1,1).type = 'text';
end

Prompt(end+1,:) = {'Reduce by PCA','doPCA'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Max number of clusters','maxclusters'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [2 9999];
Formats(end,1).size = 35;

Prompt(end+1,:) = {'Write Cluster Matrices','WriteMatrices'};
Formats(end+1,1).type = 'check';

if dosearch < 2
    Prompt(end+1,:) = {'Results-File','results'};
    Formats(end+1,1).type = 'edit';
    Formats(end,1).format = 'file';
    Formats(end,1).items = {'*.xls;*.xlsx','Select Excel-file with results'};
    Formats(end,1).limits = [0 1];
end

[settings,Cancelled] = inputsdlg(Prompt,'Kmeans clustering',Formats,settings);
if Cancelled == 1
    settings = [];
    skipprocessing = 1;
elseif isfield(settings,'searchmode')
    pause(0.2);
    if ~strcmp(settings.searchmode,'strings')
        settings.searchstring = {settings.searchmode};
    end
    switch settings.searchmode
        case 'Matrix'
            settings.searchstring = {'matrix.txt'};
    end
end