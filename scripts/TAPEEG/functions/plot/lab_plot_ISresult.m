function lab_plot_ISresult

settings = lab_set_plot_ISresult;
if isempty(settings)
    return
else
    pause(0.2);
end

if exist(settings.ris_file,'file')
    [data,header] = lab_read_data(settings.ris_file);
else
    data = [];
end
if isempty(data)
    disp('   Abort: no valid IS result file')
    return
end
if ~isfield(header,'locs') | isempty(header.locs)
    spi_file = uigetfile('*.spi','Select SPI-File with result positions');
    header.locs = lab_read_data(spi_file);
    if isempty(header.locs)
        disp('   Abort: no valid result positions file')
        return
    end
end

% Create brain mesh
if ~isempty(settings.mri_file) & exist(settings.mri_file,'file')
    cfgtmp.mrifile = settings.mri_file;
    cfgtmp.SEGcorrect = 1;
    cfgtmp.SEGprobmaps = 0;
    cfgtmp.SEGbrain = 1;
    mri = lab_segment_mri(cfgtmp);
    Mbrain = mri.gray;
    Mbrain(mri.white == 1) = 1;
    disp('    Mesh brain');
    ekernel3 = cat(3,[0 0 0; 0 1 0; 0 0 0],[0 1 0; 1 1 1; 0 1 0],[0 0 0; 0 1 0; 0 0 0]);
    tmp = spm_dilate_erode(double(Mbrain),ekernel3,'erode');
    tmp = spm_dilate_erode(double(tmp),ekernel3,'dilate');
    Mbrain(tmp==0) = 0;
    CC = bwconncomp(Mbrain);
    CC = CC.PixelIdxList';
    if size(CC,1) > 2
        for i = 3:size(CC,1)
            Mbrain(CC{i}) = 0;
        end
    end
    [Btri,Bpnt]=isosurface(Mbrain);
    [brain.tri,Bpnt] = reducepatch(Btri,Bpnt,0.1);
    Bpnt = [Bpnt(:,2) Bpnt(:,1) Bpnt(:,3)];
    if isfield(mri,'originator')
        brain.pnt = Bpnt - repmat(mri.originator,size(Bpnt,1),1) - 1;
    end
    clearvars cfgtmp mri Mbrain CC Btri Bpnt tmp
else
    brain = [];
end
if isempty(brain)
    disp('   Abort: brain mesh could not be generated')
    return
end

% plot result
plot.facecolor = [0.7 0.7 0.7];
plot.plotedges = false;
plot.plotfaces = true;
plot.alpha = 0.5;
fig1 = lab_plot_mesh(brain,plot);
set(fig1,'Name','IS result','NumberTitle','off');
P = [];

vertices = [header.locs.x' header.locs.y' header.locs.z'];
[x,y,z] = sphere(round(settings.sizedots));
x = x / 64 * settings.sizedots;
y = y / 64 * settings.sizedots;
z = z / 64 * settings.sizedots;
cmap = lab_create_cmap(settings.dotcolor);

tmp = mean(data,1);
timepoint = find(tmp == max(tmp),1);
clearvars tmp
tmp = abs(data(:,timepoint));
tmp = ceil(tmp ./ max(tmp) * 63) + 1;
for i = 1:size(vertices,1)
    P(i) = patch(surf2patch(tmp(i)*x+vertices(i,1),tmp(i)*y+vertices(i,2),tmp(i)*z+vertices(i,3)), ...
        'FaceColor',cmap(tmp(i),:),'EdgeColor','none');
end

T1 = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0,...
    'Position', [5 22 120 15],'Callback', {@(~,~)set_timepoint});
T2 = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',timepoint/size(data,2),...
    'Position', [5 5 120 15],'Callback', {@(~,~)set_timepoint});
U.Hprint = [T1 T2];
set(gcf,'UserData',U);
clearvars U
    
    function set_timepoint
        delete(P);
        numel = size(data,2);
        timepoint = round(get(T2,'Value') * numel);
        value = get(T1,'Value');
        Ttmp = abs(data(:,timepoint));
        Ttmp = ceil((Ttmp./(max(Ttmp))*(1-value)) * 63) + 1;
        Ttmp(Ttmp>64) = 64;
        for j = 1:size(vertices,1)
            P(j) = patch(surf2patch(Ttmp(j)*x+vertices(j,1),Ttmp(j)*y+vertices(j,2),Ttmp(j)*z+vertices(j,3)), ...
                'FaceColor',cmap(Ttmp(j),:),'EdgeColor','none');
        end
    end

end