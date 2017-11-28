% Plot mri using lab_plot_orthoslides
%
% written by F. Hatz 2014

function lab_plot_mri(settings,flag)

disp('Plot MRI')

if ~exist('flag','var')
    flag = false;
end
if ~exist('settings','var')
    settings = [];
end
if isfield(settings,'anatomy')
    mri = settings;
    settings = [];
elseif isnumeric(settings) & size(settings,3) > 1 & size(settings,4) == 1
    mri.anatomy = settings;
    mri.transform = eye(4);
    settings = [];
end
if ~exist('mri','var')
    if ~isfield(settings,'MRI-file')
        settings = lab_set_plot_mri(settings,flag);
        if isempty(settings)
            return
        else
            pause(0.2);
        end
    end
    mri = lab_read_mri(settings.MRI_file);
    if isempty(mri) | ~isfield(mri,'anatomy')
        return
    end
    [~,~,~,mrifileS] = lab_filename(settings.MRI_file);
else
    mrifileS = 'MRI';
end
tmp = unique(mri.anatomy(:));
if length(tmp) == 2 & max(tmp) == 1
    mri.anatomy = mri.anatomy * 255;
end
clearvars tmp

cfg = [];
if isfield(settings,'ATLAS_file') & exist(settings.ATLAS_file,'file')
    atlas = lab_match_atlas2mri(settings.ATLAS_file,settings.MRI_file);
    if ~isfield(atlas,'labels')
        tmp = setdiff(unique(atlas.anatomy),0);
        if length(tmp) == 116
            atlas.labels = lab_get_aal;
            cfg.atlas = atlas;
        else
            disp('Atlas-file disabled, no labels information')
        end
        clearvars tmp
    else
        cfg.atlas = atlas;
    end
    clearvars atlas
end
if isfield(settings,'SPI_file') & exist(settings.SPI_file,'file')
    spi = lab_read_data(settings.SPI_file);
else
    spi = [];
end
if ~isempty(spi) & isfield(spi,'x')
    if isfield(settings,'RIS_file') & exist(settings.RIS_file,'file')
        data = lab_read_data(settings.RIS_file);
    else
        data = [];
    end
    if ~isempty(data) & size(data,1) == length(spi.x)
        if size(data,2) <= 20
            Mdata = data;
        else
            sel_result = cat(1,{'Average'},cellstr(num2str((1:size(data,2))')));
            selection = listdlg('PromptString','Select Timeframe','SelectionMode','single', ...
                'ListString',sel_result,'CancelString','None','ListSize',[100 450]);
            pause(0.2);
            selection = selection - 1;
            if isempty(selection) | selection == 0
                disp('     averge data')
                Mdata = mean(data,2);
            else
                Mdata = data(:,selection);
            end
            clearvars selection sel_result
        end
    else
        Mdata = ones(length(spi.x),1);
        cfg.novalue = 1;
        if ~isfield(settings,'Interpolate') | settings.Interpolate == 0
            settings.Interpolate = 2;
        end
    end
    if ~isfield(settings,'Interpolate') | settings.Interpolate == 0
        disp('     calculate factor for interpolation (= minimal distance between solution points)')
        settings.Interpolate = round(max(min(lab_rm_diagonal(lab_distance(spi)))));
    end
    result = zeros(size(mri.anatomy,1),size(mri.anatomy,2),size(mri.anatomy,3));
    mri.result = zeros(size(mri.anatomy,1),size(mri.anatomy,2),size(mri.anatomy,3),size(Mdata,2));
    [result,T] = affine(result,mri.transform,[1 1 1]*settings.Interpolate,0,0,2);
    if isfield(mri,'originator')
        x = ceil((spi.x + mri.originator(1) + 1) / settings.Interpolate);
        y = ceil((spi.y + mri.originator(2) + 1) / settings.Interpolate);
        z = ceil((spi.z + mri.originator(3) + 1) / settings.Interpolate);
        %x = round(spi.x / settings.Interpolate) + round(mri.originator(1) / settings.Interpolate)+1;
        %y = round(spi.y / settings.Interpolate) + round(mri.originator(2) / settings.Interpolate)+1;
        %z = round(spi.z / settings.Interpolate) + round(mri.originator(3) / settings.Interpolate)+1;
    else
        x = ceil(spi.x / settings.Interpolate);
        y = ceil(spi.y / settings.Interpolate);
        z = ceil(spi.z / settings.Interpolate);
    end
    for j = 1:size(Mdata,2)
        resulttmp = result;
        for i = 1:length(x)
            resulttmp(x(i),y(i),z(i)) = Mdata(i,j);
        end
        disp(['     interpolate results ' num2str(j) ' (factor: ' num2str(settings.Interpolate) ')'])
        resulttmp = affine(resulttmp,T,[1 1 1],0,0,2);
        mri.result(:,:,:,j) = resulttmp(1:size(mri.anatomy,1),1:size(mri.anatomy,2),1:size(mri.anatomy,3));
        clearvars resulttmp
    end
    mri.spi.loc = round([spi.x(:) spi.y(:) spi.z(:)] + repmat(mri.originator,[length(spi.x) 1]) + 1);
    mri.spi.labels = spi.labels(:);
    resulttmp = mri.result(:,:,:,1);
    [~,idx] = max(resulttmp(:));
    [X,Y,Z] = ind2sub(size(resulttmp),idx(1));
    cfg.position = [X,Y,Z];
    clearvars X Y Z idx i resulttmp Mdata result x y z T
end

cfg.fhandle = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none','Name',mrifileS);
pos = get(0,'ScreenSize');
pos = [100 (pos(4)-800) 700 700];
if pos(2) < 0
    pos(2) = 100;
end
set(cfg.fhandle,'Position',pos);
lab_plot_orthoslides(mri,cfg);

end