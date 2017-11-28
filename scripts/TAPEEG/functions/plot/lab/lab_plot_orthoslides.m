% function to plot ortho-slides of mri
%
% [mri,cfg] = lab_plot_orthoslides(mri,cfg)
%
% written for TAPEEG 2014 F. Hatz

function [mri,cfg] = lab_plot_orthoslides(mri,cfg)

if ~exist('mri','var') | ~isfield(mri,'anatomy')
    mri = [];
    cfg = [];
    return
end
if ~exist('cfg','var') | ~isfield(cfg,'color')
    cfg.color = [1 0 0];
end
color = cfg.color;

Isize = max(size(mri.anatomy));
vol = zeros(Isize,Isize,Isize);
for i = 1:3
    Index1(i) = floor((Isize-size(mri.anatomy,i))/2)+1; %#ok<AGROW>
    Index2(i) = Isize - ceil((Isize-size(mri.anatomy,i))/2); %#ok<AGROW>
end
vol(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3)) = mri.anatomy;
vol = permute(vol,[2 1 3]);
volorig = vol;
if isfield(mri,'minval') & ~isempty(mri.minval)
    minval = mri.minval;
else
    minval = min(vol(:));
end
if isfield(mri,'maxval') & ~isempty(mri.maxval)
    maxval = mri.maxval;
else
    maxval = max(vol(:));
end

if isfield(mri,'result')
    maxresult = size(mri.result,4);
    numresult = 1;
    resulttmp = abs(mri.result);
    range = [min(resulttmp(:)) max(resulttmp(:))];
    result = zeros(Isize,Isize,Isize,maxresult);
    result(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3),:) = ...
        (resulttmp - min(resulttmp(:))) / (max(resulttmp(:)) - min(resulttmp(:)));
    result2 = zeros(Isize,Isize,Isize,maxresult);
    result2(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3),:) = resulttmp;
    result = permute(result,[2 1 3 4]);
else
    result = [];
end

if ~isfield(cfg,'position') | isempty(cfg.position)
    position = [round(Isize/2) round(Isize/2) round(Isize/2)];
else
    position = cfg.position;
    position = position + Index1 - 1;
    position = position([2 1 3]);
end

if isfield(cfg,'atlas') & isfield(cfg.atlas,'anatomy')
    atlas = zeros(Isize,Isize,Isize);
    atlas(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3)) = cfg.atlas.anatomy;
    atlas = permute(atlas,[2 1 3]);
    atlas_index = setdiff(unique(atlas(:)),0);
    if isfield(cfg.atlas,'labels')
        labels = cfg.atlas.labels(:);
        for i = 1:size(labels,1)
            labels{i,1} = regexprep(labels{i,1},'_',' ');
        end
    else
        for i = 1:length(atlas_index)
            labels{i,1} = ['Nr ' num2str(atlas_index(i))]; %#ok<AGROW>
        end
    end
else
    atlas = [];
end

if ~isfield(cfg,'fhandle') | ~ishandle(cfg.fhandle)
    fig = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off', ...
        'Menubar','none','Name','MRI-Ortho-Slides');
    pos = get(0,'ScreenSize');
    pos = [100 (pos(4)-800) 700 700];
    if pos(2) < 0
        pos(2) = 100;
    end
    set(fig,'Position',pos);
else
    fig = figure(cfg.fhandle);
    clf;
end
m1 = uimenu(fig,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
uimenu(m1,'Label','Close','Callback','close;');
m2 = uimenu(fig,'Label','Scroll');
m21 = uimenu(m2,'Label','Transversal');
uimenu(m21,'Label','Up','Callback',@(~,~)scroll_transversal(1),'Accelerator','1');
uimenu(m21,'Label','Down','Callback',@(~,~)scroll_transversal(-1),'Accelerator','2');
m22 = uimenu(m2,'Label','Sagittal');
uimenu(m22,'Label','Up','Callback',@(~,~)scroll_saggital(1),'Accelerator','3');
uimenu(m22,'Label','Down','Callback',@(~,~)scroll_sagittal(-1),'Accelerator','4');
m23 = uimenu(m2,'Label','Coronar');
uimenu(m23,'Label','Up','Callback',@(~,~)scroll_coronar(1),'Accelerator','5');
uimenu(m23,'Label','Down','Callback',@(~,~)scroll_coronar(-1),'Accelerator','6');
if isempty(result)
    m3 = uimenu(fig,'Label','Edit');
    uimenu(m3,'Label','Goto Position','Callback',@(~,~)goto_position);
    uimenu(m3,'Label','Set Min/Max','Callback',@(~,~)set_min_max);
    uimenu(m3,'Label','Anisotropic Filtering','Callback',@(~,~)filter_mri);
    uimenu(m3,'Label','Correct orientation to RAS','Callback',@(~,~)correct_orient);
    uimenu(m3,'Label','Set Landmarks','Callback',@(~,~)set_landmarks);
    m4 = uimenu(m3,'Label','Originator');
    uimenu(m4,'Label','Set','Callback',@(~,~)set_originator);
    uimenu(m4,'Label','Goto','Callback',@(~,~)goto_originator);
    uimenu(m1,'Label','Save MRI','Callback',@(~,~)save_mri,'Separator','on');
else
    m3 = uimenu(fig,'Label','Select');
    uimenu(m3,'Label','Maximal Result','Callback',@(~,~)goto_maxresult);
    if isfield(mri,'spi') & isfield(mri.spi,'loc')
        if ~isfield(mri.spi,'labels')
            mri.spi.labels = cellstr(num3str((1:size(mri.spi.loc,1))'));
        end
        uimenu(m3,'Label','Select Position','Callback',@(~,~)select_position);
    end
    if size(result,4) > 1
        uimenu(m3,'Label','Select Result','Callback',@(~,~)select_result);
    end
end

IMx = calc_ximage(position(1));
IMy = calc_yimage(position(2));
IMz = calc_zimage(position(3));

U = get(gcf,'UserData');
if ~isfield(U,'Hprint')
    U.Hprint = [];
end

Az = axes('Position',[0.015 0.535 0.45 0.45],'Color',[0 0 0]);
Hz = image(IMz,'ButtonDownFcn',@(~,~)refresh_image_z);
set(gca,'Visible','off','PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1],'Color',[0 0 0], ...
    'YDir','normal','XDir','reverse');
hold on
Pz = plot_crosshair(position(2),position(1),@(~,~)refresh_image_z);
text(0,0,'\leftarrow X ','FontSize',9,'FontUnits','normalized','HorizontalAlignment', ...
    'right','VerticalAlignment','top')
text(0,0,' Y \rightarrow','FontSize',9,'FontUnits','normalized','HorizontalAlignment', ...
    'left','VerticalAlignment','top','Rotation',90)

Ax = axes('Position',[0.015 0.04 0.45 0.45],'Color',[0 0 0]);
Hx = image(IMx,'ButtonDownFcn',@(~,~)refresh_image_x);
set(gca,'Visible','off','PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1],'Color',[0 0 0], ...
    'YDir','normal','XDir','reverse');
hold on
Px = plot_crosshair(position(2),position(3),@(~,~)refresh_image_x);
text(0,0,'\leftarrow X ','FontSize',9,'FontUnits','normalized','HorizontalAlignment', ...
    'right','VerticalAlignment','top')
text(0,0,' Z \rightarrow','FontSize',9,'FontUnits','normalized','HorizontalAlignment', ...
    'left','VerticalAlignment','top','Rotation',90)

Ay = axes('Position',[0.51 0.535 0.45 0.45],'Color',[0 0 0]);
Hy = image(IMy,'ButtonDownFcn',@(~,~)refresh_image_y);
set(gca,'Visible','off','PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1],'Color',[0 0 0], ...
    'YDir','normal','XDir','reverse');
hold on
Py = plot_crosshair(position(1),position(3),@(~,~)refresh_image_y);
text(0,0,'\leftarrow Y ','FontSize',9,'FontUnits','normalized','HorizontalAlignment', ...
    'right','VerticalAlignment','top')
text(0,0,' Z \rightarrow','FontSize',9,'FontUnits','normalized','HorizontalAlignment', ...
    'left','VerticalAlignment','top','Rotation',90)

axes('Position',[0.525 0.025 0.45 0.45],'Color',[1 1 1],'Visible','off');
pos = position;
pos = pos([2 1 3]);
pos = pos - Index1 + 1;
if isfield(mri,'originator')
    pos = pos - mri.originator - 1;
end
T1 = text(0.1,0.9,['Position:  ' num2str(pos(1)) '  ' num2str(pos(2)) '  ' num2str(pos(3))], ...
    'FontSize',9,'FontUnits','normalized','UserData',position);
if ~isempty(atlas)
    Iatlas = atlas(position(2),position(1),position(3));
    if Iatlas ~= 0
        Iatlas = find(atlas_index==Iatlas);
        if ~isempty(Iatlas)
            T2 = text(0.1,0.7,['Label:  ' labels{Iatlas}],'FontSize',9,'FontUnits','normalized');
        else
            T2 = text(0.1,0.7,'Label:  ','FontSize',9,'FontUnits','normalized');
        end
    else
        T2 = text(0.1,0.7,'Label:  ','FontSize',9,'FontUnits','normalized');
    end
else
    T2 = [];
end

if ~isempty(result) & ~isfield(cfg,'novalue')
    Rvalue = result2(position(2),position(1),position(3),numresult);
    T3 = text(0.1,0.5,['Value:  ' num2str(Rvalue)],'FontSize',9,'FontUnits','normalized');
    text(0.1,0.42,['Min:  ' num2str(range(1))],'FontSize',9,'FontUnits','normalized');
    text(0.1,0.34,['Max:  ' num2str(range(2))],'FontSize',9,'FontUnits','normalized');
else
    T3 = [];
end

SliderZ = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0.5,...
    'Units','normalized','Position', [0.14 0.51 0.2 0.025],'Callback', {@slider_z});
SliderX = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0.5,...
    'Units','normalized','Position', [0.14 0.015 0.2 0.025],'Callback', {@slider_x});
SliderY = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0.5,...
    'Units','normalized','Position', [0.635 0.51 0.2 0.025],'Callback', {@slider_y});
U.Hprint = [U.Hprint SliderX SliderY SliderZ];
set(gcf,'UserData',U);

if isfield(cfg,'showorient') & cfg.showorient == true
    correct_orient;
end

if isfield(cfg,'setlandmarks') & cfg.setlandmarks == true
    set_landmarks;
end

if isfield(cfg,'setminmax') & cfg.setminmax == true
    set_min_max;
end

% inline scripts
    function refresh_image_z
        Pos = get(Az,'CurrentPoint');
        position(2) = round(Pos(1,1));
        position(1) = round(Pos(1,2));
        position(position<1) = 1;
        position(position>Isize) = Isize;
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        refresh_image;
    end

    function refresh_image_x
        Pos = get(Ax,'CurrentPoint');
        position(2) = round(Pos(1,1));
        position(3) = round(Pos(1,2));
        position(position<1) = 1;
        position(position>Isize) = Isize;
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function refresh_image_y
        Pos = get(Ay,'CurrentPoint');
        position(1) = round(Pos(1,1));
        position(3) = round(Pos(1,2));
        position(position<1) = 1;
        position(position>Isize) = Isize;
        IMx = calc_ximage(position(1));
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function scroll_transversal(skip)
        position(3) = position(3) + skip;
        position(position<1) = 1;
        position(position>Isize) = Isize;
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function slider_x(hObj,~)
        position(1) = round(get(hObj,'Value') * (Isize-1))+1;
        IMx = calc_ximage(position(1));
        refresh_image;
    end

    function slider_y(hObj,~)
        position(2) = round(get(hObj,'Value') * (Isize-1))+1;
        IMy = calc_yimage(position(2));
        refresh_image;
    end

    function slider_z(hObj,~)
        position(3) = round(get(hObj,'Value') * (Isize-1))+1;
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function scroll_saggital(skip)
        position(2) = position(2) + skip;
        position(position<1) = 1;
        position(position>Isize) = Isize;
        IMy = calc_yimage(position(2));
        refresh_image;
    end

    function scroll_coronar(skip)
        position(1) = position(1) + skip;
        position(position<1) = 1;
        position(position>Isize) = Isize;
        IMx = calc_ximage(position(1));
        refresh_image;
    end

    function refresh_image
        set(Hx,'CData',IMx);
        set(Hy,'CData',IMy);
        set(Hz,'CData',IMz);
        set(Px(2),'XData',[position(2) position(2)]);
        set(Px(1),'YData',[position(3) position(3)]);
        set(Py(2),'XData',[position(1) position(1)]);
        set(Py(1),'YData',[position(3) position(3)]);
        set(Pz(1),'YData',[position(1) position(1)]);
        set(Pz(2),'XData',[position(2) position(2)]);
        set(SliderX,'Value',(position(1)-1)/Isize);
        set(SliderY,'Value',(position(2)-1)/Isize);
        set(SliderZ,'Value',(position(3)-1)/Isize);
        pos2 = position ([2 1 3]);
        pos2 = pos2 - Index1 + 1;
        if isfield(mri,'originator')
            pos2 = pos2 - mri.originator - 1;
        end
        set(T1,'String',['Position:  ' num2str(pos2(1)) '  ' num2str(pos2(2)) '  ' num2str(pos2(3))],'UserData',position);
        if ~isempty(T2)
            Iatlas = atlas(position(1),position(2),position(3));
            if Iatlas ~= 0
                Iatlas = find(atlas_index==Iatlas);
                if ~isempty(Iatlas)
                    set(T2,'String',['Label:  ' labels{Iatlas}]);
                else
                    set(T2,'String','Label:  ');
                end
            else
                set(T2,'String','Label:  ');
            end
        end
        if ~isempty(T3)
            Rvalue = result2(position(2),position(1),position(3),numresult);
            set(T3,'String',['Value:  ' num2str(Rvalue)]);
        end
    end

    function IMx = calc_ximage(xpos)
        xpos = round(xpos);
        if xpos > size(vol,1)
            xpos = size(vol,1);
        elseif xpos < 1
            xpos = 1;
        end
        voltmp = (vol - minval) / ((maxval-minval)/255);
        Vtmp = repmat(permute(voltmp(xpos,:,:),[3 2 1]),[1 1 3]);
        Vtmp(Vtmp<0) = 0;
        Vtmp(Vtmp>255) = 255;
        if ~isempty(result)
            Vtmp2 = repmat(permute(result(xpos,:,:,numresult),[3 2 1]),[1 1 3]);
            Vtmp3(:,:,1) = Vtmp2(:,:,1) * color(1) * 255;
            Vtmp3(:,:,2) = Vtmp2(:,:,2) * color(2) * 255;
            Vtmp3(:,:,3) = Vtmp2(:,:,3) * color(3) * 255;
            IMx = uint8(Vtmp .* (1-Vtmp2) + Vtmp3 .* Vtmp2);
        else
            IMx = uint8(Vtmp);
        end
    end

    function IMy = calc_yimage(ypos)
        ypos = round(ypos);
        if ypos > size(vol,2)
            ypos = size(vol,2);
        elseif ypos < 1
            ypos = 1;
        end
        voltmp = (vol - minval) / ((maxval-minval)/255);
        Vtmp = repmat(permute(voltmp(:,ypos,:),[3 1 2]),[1 1 3]);
        Vtmp(Vtmp<0) = 0;
        Vtmp(Vtmp>255) = 255;
        if ~isempty(result)
            Vtmp2 = repmat(permute(result(:,ypos,:,numresult),[3 1 2]),[1 1 3]);
            Vtmp3(:,:,1) = Vtmp2(:,:,1) * color(1) * 255;
            Vtmp3(:,:,2) = Vtmp2(:,:,2) * color(2) * 255;
            Vtmp3(:,:,3) = Vtmp2(:,:,3) * color(3) * 255;
            IMy = uint8(Vtmp .* (1-Vtmp2) + Vtmp3 .* Vtmp2);
        else
            IMy = uint8(Vtmp);
        end
    end

    function IMz = calc_zimage(zpos)
        zpos = round(zpos);
        if zpos > size(vol,3)
            zpos = size(vol,3);
        elseif zpos < 1
            zpos = 1;
        end
        voltmp = (vol - minval) / ((maxval-minval)/255);
        Vtmp = repmat(voltmp(:,:,zpos),[1 1 3]);
        Vtmp(Vtmp<0) = 0;
        Vtmp(Vtmp>255) = 255;
        if ~isempty(result)
            Vtmp2 = repmat(result(:,:,zpos,numresult),[1 1 3]);
            Vtmp3(:,:,1) = Vtmp2(:,:,1) * color(1) * 255;
            Vtmp3(:,:,2) = Vtmp2(:,:,2) * color(2) * 255;
            Vtmp3(:,:,3) = Vtmp2(:,:,3) * color(3) * 255;
            IMz = uint8(Vtmp .* (1-Vtmp2) + Vtmp3 .* Vtmp2);
        else
            IMz = uint8(Vtmp);
        end
    end

    function P = plot_crosshair(y,x,fc)
        Ysize = get(gca,'YLim');
        Xsize = get(gca,'XLim');
        P(1) = plot(Ysize,[x x],'-b','ButtonDownFcn',fc);
        P(2) = plot([y y],Xsize,'-b','ButtonDownFcn',fc);
    end

    function filter_mri
        settings = lab_set_mrifiltering;
        pause(0.2);
        if ~isempty(settings)
            voltmp = vol;
            voltmp(voltmp>maxval) = maxval;
            voltmp = voltmp - minval;
            voltmp(voltmp<0) = 0;
            mri.anatomy = permute(voltmp(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3)),[2 1 3]);
            mri = lab_filter_mri(mri,settings);
            vol = zeros(Isize,Isize,Isize);
            vol(Index1(2):Index2(2),Index1(1):Index2(1),Index1(3):Index2(3)) = mri.anatomy;
            vol = permute(vol,[2 1 3]);
            volorig = vol;
            IMx = calc_ximage(position(1));
            IMy = calc_yimage(position(2));
            IMz = calc_zimage(position(3));
            refresh_image;
        end
    end

    function save_mri
         [Filename,Filepath] = uiputfile('*.hdr','Store MRI-file');
         mri.anatomy = permute(vol(Index1(2):Index2(2),Index1(1):Index2(1),Index1(3):Index2(3)),[2 1 3]);
         lab_write_hdr(fullfile(Filepath,Filename),mri);
    end

    function select_position
        selection = listdlg('PromptString','Select Position','SelectionMode','single', ...
            'ListString',mri.spi.labels,'CancelString','None','ListSize',[100 450]);
        if ~isempty(selection)
            position = mri.spi.loc(selection,:);
            position = position + Index1 - 1;
            position = position([2 1 3]);
            IMx = calc_ximage(position(1));
            IMy = calc_yimage(position(2));
            IMz = calc_zimage(position(3));
            refresh_image;
        end
    end

    function select_result
        sel_result = cellstr(num2str((1:maxresult)'));
        selection = listdlg('PromptString','Select Result','SelectionMode','single', ...
            'ListString',sel_result,'CancelString','None','ListSize',[100 450]);
        if ~isempty(selection)
            numresult = selection;
            IMx = calc_ximage(position(1));
            IMy = calc_yimage(position(2));
            IMz = calc_zimage(position(3));
            refresh_image;
        end
    end

    function goto_maxresult
        Tresult = result(:,:,:,numresult);
        [~,idx] = max(Tresult(:));
        [X,Y,Z] = ind2sub(size(Tresult),idx(1));
        position = [X,Y,Z];
        clearvars X Y Z idx Tresult
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function correct_orient
        [mri,cfg.orient] = lab_correct_mri_orient(mri);
        Isize = max(size(mri.anatomy));
        vol = zeros(Isize,Isize,Isize);
        for j = 1:3
            Index1(j) = floor((Isize-size(mri.anatomy,j))/2)+1;
            Index2(j) = Isize - ceil((Isize-size(mri.anatomy,j))/2);
        end
        vol(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3)) = mri.anatomy;
        vol = permute(vol,[2 1 3]);
        volorig = vol;
        atlas = [];
        if ~isempty(T2)
            delete(T2)
            T2 = [];
        end
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end
    
    function goto_position
        settings.pos = position;
        settings.pos = settings.pos([2 1 3]);
        settings.pos = settings.pos - Index1 + 1;
        if isfield(mri,'originator')
            settings.pos = settings.pos - mri.originator - 1;
        end
        Prompt = {'Position','pos'};
        Formats.type = 'edit';
        Formats.format = 'vector';
        Formats.limits = [-inf inf];
        Formats.items = [];
        [settings,Cancelled] = inputsdlg(Prompt,'Position',Formats,settings);
        if Cancelled ~= 1
            if size(settings.pos,2) ~= 3
                return
            end
            settings.pos = round(settings.pos);
            if isfield(mri,'originator')
                settings.pos = settings.pos + mri.originator + 1;
            end
            settings.pos = settings.pos + Index1 - 1;
            settings.pos = settings.pos([2 1 3]);
            if min(settings.pos) < 1 | max(settings.pos) > Isize
                return
            end
            position = settings.pos;
            IMx = calc_ximage(position(1));
            IMy = calc_yimage(position(2));
            IMz = calc_zimage(position(3));
            refresh_image;
        end
    end

    function set_landmarks
        if ~isfield(mri,'landmarks') | isempty(mri.landmarks)
            answer = questdlg('Standard or Define by electrodes-file','Select landmarks','File','Standard','Standard');
            selection = {};
            if strcmp(answer,'File')
                [tmp,~,LOCS] = lab_plot_locs([],1);
                if ~isempty(tmp)
                    selection = LOCS.labels(1,tmp);
                end
            end
            if isempty(selection)
                selection = {'Nasion','Inion','LPA','RPA','Cz'};
            end
            for S = 1:length(selection)
                landmarks(S).name = selection{S}; %#ok<AGROW>
                landmarks(S).pnt = []; %#ok<AGROW>
                landmarks(S).mode = 'fixed'; %#ok<AGROW>
                landmarks(S).calcfactor = []; %#ok<AGROW>
                landmarks(S).calc1 = []; %#ok<AGROW>
                landmarks(S).calc2 = []; %#ok<AGROW>
            end
        else
            landmarks = mri.landmarks;
        end
        Prompt = cell(0,2);
        Formats = [];
        list = {' '};
        for j = 1:length(landmarks)
            list{end+1,1} = landmarks(j).name; %#ok<AGROW>
        end
        for j = 1:length(landmarks)
            Prompt(end+1,:) = {landmarks(j).name,''}; %#ok<AGROW>
            Formats(end+1,1).type = 'text'; %#ok<AGROW>
            
            settings.(['pnt' num2str(j)]) = landmarks(j).pnt;
            Prompt(end+1,:) = {'',['pnt' num2str(j)]}; %#ok<AGROW>
            Formats(end+1,1).type = 'edit'; %#ok<AGROW>
            Formats(end,1).format = 'vector';
            Formats(end,1).limits = [-inf inf];
            Formats(end,1).size = 70;
            
            Prompt(end+1,:) = {'Set',''}; %#ok<AGROW>
            Formats(end+1,1).type = 'button'; %#ok<AGROW>
            Formats(end,1).style = 'pushbutton';
            Formats(end,1).size = [50 25];
            Formats(end,1).callback = {@set_pnt,['pnt' num2str(j)],['pnt' num2str(j)]};
            
            Prompt(end+1,:) = {'Goto',''}; %#ok<AGROW>
            Formats(end+1,1).type = 'button'; %#ok<AGROW>
            Formats(end,1).style = 'pushbutton';
            Formats(end,1).size = [50 25];
            Formats(end,1).callback = {@goto_pnt,[],['pnt' num2str(j)]};
                        
            settings.(['mode' num2str(j)]) = landmarks(j).mode;
            Prompt(end+1,:) = {'Mode',['mode' num2str(j)]}; %#ok<AGROW>
            Formats(end+1,1).type = 'list'; %#ok<AGROW>
            Formats(end,1).style = 'popupmenu';
            Formats(end,1).format = 'input';
            Formats(end,1).items = {'fixed','correct'};
            Formats(end,1).callback = {@set_mode,{['calcfactor' num2str(j)], ...
                ['calc1_' num2str(j)],['calc2_' num2str(j)]},['mode' num2str(j)], ...
                ['calcfactor' num2str(j)],['calc1_' num2str(j)],['calc2_' num2str(j)]};
            
            settings.(['calcfactor' num2str(j)]) = landmarks(j).calcfactor;
            Prompt(end+1,:) = {'Percent',['calcfactor' num2str(j)]}; %#ok<AGROW>
            Formats(end+1,1).type = 'edit'; %#ok<AGROW>
            Formats(end,1).format = 'integer';
            Formats(end,1).limits = [1 99];
            Formats(end,1).size = 30;
            
            settings.(['calc1_' num2str(j)]) = landmarks(j).calc1;
            Prompt(end+1,:) = {'from',['calc1_' num2str(j)]}; %#ok<AGROW>
            Formats(end+1,1).type = 'list'; %#ok<AGROW>
            Formats(end,1).style = 'popupmenu';
            Formats(end,1).format = 'input';
            Formats(end,1).items = list;
            
            settings.(['calc2_' num2str(j)]) = landmarks(j).calc2;
            Prompt(end+1,:) = {'to',['calc2_' num2str(j)]}; %#ok<AGROW>
            Formats(end+1,1).type = 'list'; %#ok<AGROW>
            Formats(end,1).style = 'popupmenu';
            Formats(end,1).format = 'input';
            Formats(end,1).items = list;
        end
        Options.WindowStyle = 'normal';
        Options.ForceWindowStyle = true;
        Options.Dim = 8;
        [settings,Cancelled] = inputsdlg(Prompt,'Set Landmarks',Formats,settings,Options);
        if Cancelled == 0
            for j = 1:length(landmarks)
                landmarks(j).pnt = settings.(['pnt' num2str(j)]); %#ok<AGROW>
                landmarks(j).mode = settings.(['mode' num2str(j)]); %#ok<AGROW>
                landmarks(j).calcfactor = settings.(['calcfactor' num2str(j)]); %#ok<AGROW>
                landmarks(j).calc1 = settings.(['calc1_' num2str(j)]); %#ok<AGROW>
                landmarks(j).calc2 = settings.(['calc2_' num2str(j)]); %#ok<AGROW>
            end
            mri.landmarks = landmarks;
            mri.storelandmarks = true;
        else
            mri.storelandmarks = false;
        end
    end
    
    function set_originator
        mri.originator = position([2 1 3]) - Index1;
        if isfield(mri,'transform')
            mri.transform(1:3,4) = -mri.originator;
        end
        refresh_image;
    end
    
    function goto_originator
        position = round(mri.originator + Index1);
        position = position([2 1 3]);
        if ~isempty(result)
            result(position(1),position(2),position(3),numresult) = 1;
        end
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end
    
    function pnt = set_pnt(pnt)
        if isempty(result)
            result = zeros(Isize,Isize,Isize,1);
            maxresult = 1;
            numresult = 1;
            result2 = result;
        end
        if ~isempty(pnt)
            pnt2 = round(pnt);
            if isfield(mri,'originator')
                pnt2 = pnt2 + mri.originator + 1;
            end
            pnt2 = pnt2 + Index1 - 1;
            pnt2 = pnt2([2 1 3]);
            pnt2(pnt2<1) = 1;
            pnt2(pnt2>Isize) = Isize;
            result(pnt2(1),pnt2(2),pnt2(3),numresult) = 0;
        end
        pnt = position([2 1 3]);
        pnt = pnt - Index1 + 1;
        if isfield(mri,'originator')
            pnt = pnt - mri.originator - 1;
        end
        position = round(position);
        result(position(1),position(2),position(3),numresult) = 1;
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function goto_pnt(pnt)
        if isfield(mri,'originator')
            pnt = pnt + mri.originator + 1;
        end
        pnt = pnt + Index1 - 1;
        position = pnt([2 1 3]);
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end

    function set_min_max
        settings.minval = minval;
        settings.maxval = maxval;
        settings.excludezeros = true;
        
        Prompt = cell(0,2);
        Formats = [];
        
        Prompt(end+1,:) = {'Minimal value','minval'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'float';
        Formats(end,1).limits = [0 inf];
        Formats(end,1).size = 40;
        
        Prompt(end+1,:) = {'Set',''};
        Formats(end+1,1).type = 'button';
        Formats(end,1).style = 'pushbutton';
        Formats(end,1).size = [40 25];
        Formats(end,1).callback = {@set_contrast,[],'@ALL'};
        
        Prompt(end+1,:) = {'Maximal value','maxval'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'float';
        Formats(end,1).limits = [0 inf];
        Formats(end,1).size = 40;
        
        Prompt(end+1,:) = {'Set',''};
        Formats(end+1,1).type = 'button';
        Formats(end,1).style = 'pushbutton';
        Formats(end,1).size = [40 25];
        Formats(end,1).callback = {@set_contrast,[],'@ALL'};
        
        Formats(end+1,1).type = 'none';
        Formats(end,1).span = [1 2];
        
        Prompt(end+1,:) = {'Histogram',''};
        Formats(end+1,1).type = 'button';
        Formats(end,1).style = 'pushbutton';
        Formats(end,1).size = [100 25];
        Formats(end,1).callback = {@lab_plot_histogram,'',volorig,'@ALL'};
        
        Prompt(end+1,:) = {'Exclude zeros','excludezeros'};
        Formats(end+1,1).type = 'check';
        
        Options.WindowStyle = 'normal';
        Options.ForceWindowStyle = true;
        [settings,Cancelled] = inputsdlg(Prompt,'Set Min/Max Values',Formats,settings,Options);
        if Cancelled == 0
            set_contrast(settings)
        end
    end
    function set_contrast(settings)
        minval = settings.minval;
        maxval = settings.maxval;
        voltmp = volorig;
        voltmp(voltmp>maxval) = maxval;
        voltmp = voltmp - minval;
        voltmp(voltmp<0) = 0;
        vol = voltmp;
        mri.anatomy = permute(vol(Index1(1):Index2(1),Index1(2):Index2(2),Index1(3):Index2(3)),[2 1 3]);
        mri.minval = minval;
        mri.maxval = maxval;
        IMx = calc_ximage(position(1));
        IMy = calc_yimage(position(2));
        IMz = calc_zimage(position(3));
        refresh_image;
    end
end

function [factor,calc1,calc2] = set_mode(mode,factor,calc1,calc2)
  if strcmp(mode,'fixed')
      factor = [];
      calc1 = ' ';
      calc2 = ' ';
  elseif isempty(factor)
      factor = 50;
  end
end

