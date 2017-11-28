% Plot lectrodes locations on a scalp-mesh (optionally + brain-mesh)
% Electrodes are plotted as small spheres
% Mappings-files (*.xls) can be loaded and modified
%
% [Indexed,Mappings] = lab_plot_locs3d(enablesave)
%
% written by F. Hatz

function [Indexed,Mappings] = lab_plot_locs3d(settings,enablesave)

Mappings = [];
Indexed = [];

if ~exist('enablesave','var')
    enablesave = 0;
end
if ~exist('settings','var')
    settings = [];
end

if isempty(settings) | ~isfield(settings,'loc_file')
    [settings,skipprocessing] = lab_set_plot_locs3d(settings);
    if skipprocessing == 1
        return
    end
end

% Segment and mesh mrifile
if ~isempty(settings.mrifile) & exist(settings.mrifile,'file')
    cfgtmp.mrifile = settings.mrifile;
    cfgtmp.SEGcorrect = 1;
    cfgtmp.SEGprobmaps = 0;
    cfgtmp.SEGbrain = 1;
    mri = lab_segment_mri(cfgtmp);
    mri.brain = mri.gray;
    mri.brain(mri.white == 1) = 1;
    cfgtmp = [];
    cfgtmp.nvert = 5000;
    [bnd,cfgtmp] = lab_mesh_segm(mri,cfgtmp);
    if isfield(cfgtmp,'tissue')
        tmp = find(strcmp(cfgtmp.tissue,'scalp'));
        if isempty(tmp)
            scalp = bnd(1);
        else
            scalp = bnd(tmp);
        end
        tmp = find(strcmp(cfgtmp.tissue,'brain'));
        if isempty(tmp)
            brain = bnd(end);
        else
            brain = bnd(tmp);
        end
    else
        scalp = bnd(1);
        brain = bnd(end);
    end
else
    scalp = [];
    brain = [];
end

% read locs
if ischar(settings.locfile) & exist(settings.locfile,'file')
    locs = lab_read_locs(settings.locfile);
elseif isstruct(settings.locfile) & isfield(settings.locfile,'x')
    locs = settings.loc_file;
else
    locs = [];
end
if isempty(locs)
    return
end

% create figure
fig1 = figure('Color',[1 1 1],'Menubar','none');
axes1 = axes;
set(axes1,'Visible','off','Parent',fig1,...
    'Position',[.01 .02 .99 .95],'xtick',[],'ytick',[], ...
    'PlotBoxAspectRatio',[1 1 1],'DataAspectRatio',[1 1 1]);
view(axes1,[-90 0]); % left
Hlight(1) = light('Parent',axes1,'Position',[0 0 1]);
Hlight(2) = light('Parent',axes1,'Position',[0 0 -1]);
U.figuretype = 'plotmesh';
set(fig1,'Userdata',U)
clearvars U

% Create menu
m1 = uimenu(fig1,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
if enablesave == 1
    uimenu(m1,'Label','Save','Callback',@(~,~)do_save,'Separator','on');
    uimenu(m1,'Label','Cancel','Callback',@(~,~)do_cancel);
else
    uimenu(m1,'Label','Close','Callback','close;');
end
m2 = uimenu(fig1,'Label','Edit');
uimenu(m2,'Label','Add Mesh','Callback',@(~,~)lab_plot_mesh);
if ~isempty(scalp)
    uimenu(m2,'Label','Edit scalp','Callback',@(~,~)edit_scalp);
    if settings.plotbrain == true
        uimenu(m2,'Label','Edit brain','Callback',@(~,~)edit_brain);
    end
end
uimenu(m2,'Label','Black/White Background','Callback',@(~,~)set_background);
uimenu(m2,'Label','Rotate on/off','Callback',@(~,~)switch_rotate,'Accelerator','1');
m3 = uimenu(fig1,'Label','View');
uimenu(m3,'Label','Top','Callback',@(~,~)view_top);
uimenu(m3,'Label','Bottom','Callback',@(~,~)view_bottom);
uimenu(m3,'Label','Left','Callback',@(~,~)view_left);
uimenu(m3,'Label','Right','Callback',@(~,~)view_right);
uimenu(m3,'Label','Front','Callback',@(~,~)view_front);
uimenu(m3,'Label','Back','Callback',@(~,~)view_back);
m4 = uimenu(fig1,'Label','Mapping');
uimenu(m4,'Label','Select','Callback',@(~,~)select_mapping);
uimenu(m4,'Label','Edit','Callback',@(~,~)edit_mapping);
uimenu(m4,'Label','Add','Callback',@(~,~)add_mapping);
uimenu(m4,'Label','Delete','Callback',@(~,~)del_mapping);
uimenu(m4,'Label','Legend on/off','Callback',@(~,~)show_mappinglegend,'Separator','on');
uimenu(m4,'Label','Load','Callback',@(~,~)import_mapping,'Separator','on');
uimenu(m4,'Label','Save','Callback',@(~,~)save_mapping);

% plot result using lab_plot_mesh
plot.fhandle = fig1;
plot.facecolor = settings.facecolor;
if ~isempty(scalp)
    plot.plotedges = false;
    plot.plotfaces = true;
    plot.plotdots = false;
    plot.alpha = settings.alpha;
    [~,Hscalp] = lab_plot_mesh(scalp,plot);
    set(fig1,'Name','Locs 3D','NumberTitle','off');
    
    if settings.plotbrain == true
        plot.facecolor = settings.Bfacecolor;
        plot.alpha = settings.Balpha;
        [~,Hbrain] = lab_plot_mesh(brain,plot);
    else
        Hbrain = [];
    end
else
    Hscalp = [];
    Hbrain = [];
end

% Plot Locs
color = settings.dotscolor;
plot.alpha = 1;
plot.sizedots = settings.sizedots;
plot.plotlabels = true;
plot.labelsize = settings.labelsize;
plot.labeldistance = settings.labeldistance;
plot.plotdots = true;
plot.dotscolor = color;
[~,Hlocs,Htext] = lab_plot_mesh(locs,plot);
for i = 1:size(Hlocs,2)
    set(Hlocs(i),'UserData',i,'ButtonDownFcn',@set_channel)
    set(Htext(i),'UserData',i,'ButtonDownFcn',@set_channel)
end

U.Hprint = [];
Hrot = uicontrol('style', 'checkbox', 'units', 'pixels','position', [2,2,80,25],...
    'string','Rotate','fontsize',10,'Callback','rotate3d','Value',true,...
    'TooltipString','Rotate on/off','BackgroundColor',[1 1 1]);
U.Hprint(1,end+1) = Hrot;
set(fig1,'UserData',U);
if enablesave == 1
    U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'normalized',...
        'position', [0.875 0.01 0.12 0.05],...
        'string','Save','fontsize',12,'Callback',@(~,~)do_save,...
        'TooltipString','Save');
    U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'normalized',...
        'position', [0.75 0.01 0.12 0.05],...
        'string','Cancel','fontsize',12,'Callback',@(~,~)do_cancel,...
        'TooltipString','Cancel');
    set(Fig,'UserData',U);
    uiwait
end

numchans = size(locs.x,2);
Nmapping = 1;
mappings(1).Name = 'Selection';
mappings(1).Color = [1 0 0];
indexed = zeros(1,numchans);
Indexed = [];
Hlegend = [];

    function set_channel(H,H2)
        Nchan = get(H,'UserData');
        T = get(Hlocs(Nchan));
        if min(T.FaceColor - color == 0) == 1
            indexed(Nmapping,Nchan) = 1;
            set(Hlocs(Nchan),'FaceColor',mappings(1,Nmapping).Color);
        else
            indexed(:,Nchan) = 0;
            set(Hlocs(Nchan),'FaceColor',color);
        end
    end
    function do_save
        if size(indexed,1) > 1
            for I = 1:size(indexed,1)
                Indexed{I} = find(indexed(I,:)==1);
            end
        else
            Indexed = find(indexed==1);
        end
        Mappings = mappings;
        close;
    end
    function do_cancel
        Indexed = [];
        Mappings = [];
        close;
    end
    function switch_rotate
        T = get(Hrot,'Value');
        if T == 1
            set(Hrot,'Value',0);
            rotate3d off
        else
            set(Hrot,'Value',1);
            rotate3d on
        end
    end
    function add_mapping
        if isfield(mappings,'active')
            return
        end
        Nmapping = size(mappings,2)+1;
        indexed(Nmapping,:) = zeros(1,numchans);
        mappings(1,Nmapping).Name = ['Mapping' num2str(Nmapping)];
        Ctmp2 = get_colormap(Nmapping);
        for j = 1:size(indexed,1)
            mappings(1,j).Color = Ctmp2(j,:);
            for M = find(indexed(j,:)==1)
                set(Hlocs(M),'FaceColor',mappings(1,j).Color);
            end
        end
        clearvars Ctmp2;
        mappings(1,Nmapping) = inputsdlg(mappings(1,Nmapping),'Mappings');
    end
    function edit_mapping
        mappings2 = mappings;
        items{3}.type = 'button';
        items{3}.style = 'pushbutton';
        items{3}.size = [80 25];
        items{3}.callback = {@set_mapping,'@ALL','@ALL','Select'};
        for j = 1:length(mappings2)
            mappings2(1,j).Select = j;
            if j == Nmapping
                mappings2(1,j).Active = true;
            else
                mappings2(1,j).Active = false;
            end
        end
        Options.WindowStyle = 'normal';
        Options.CancelButton = 'off';
        [mappings2,Cancelled] = inputsdlg(mappings2,'Mappings',items,[],Options);
        if Cancelled == 0
            mappings = mappings2;
        end
        for j = 1:size(indexed,1)
            for M = find(indexed(j,:)==1)
                set(Hlocs(M),'FaceColor',mappings(1,j).Color);
            end
        end
    end
    function settings = set_mapping(settings,select)
        Nmapping = select;
        Fvar = fieldnames(settings);
        for j = 1:length(Fvar)
            if ~isempty(strfind(Fvar{j},'Active'))
                Ftmp = str2num(Fvar{j}(2:5)); %#ok<ST2NM>
                if Ftmp == Nmapping
                    settings.(Fvar{j}) = true;
                else
                    settings.(Fvar{j}) = false;
                end
            end
        end
    end
    function select_mapping
        if isfield(mappings,'active')
            return
        end
        for j = 1:size(mappings,2)
            list{j} = mappings(1,j).Name; %#ok<AGROW>
        end
        selection = listdlg('PromptString','Select','Name','Mapping','SelectionMode','single', ...
            'ListString',list,'ListSize',[250 280]);
        if ~isempty(selection)
            Nmapping = selection;
        end
    end
    function del_mapping
        if isfield(mappings,'active')
            return
        end
        if size(mappings,2) == 1
            mappings.Name = 'Mapping1';
            indexed = zeros(1,numchans);
            for j = 1:size(Lhandle,2)
                set(Hlocs(j),'FaceColor',color);
            end
        else
            include = setdiff(1:size(mappings,2),Nmapping);
            for j = find(indexed(Nmapping,:)==1)
                set(Hlocs(j),'FaceColor',color);
            end
            mappings = mappings(1,include);
            indexed = indexed(include,:);
            Ctmp2 = get_colormap(size(mappings,2));
            for j = 1:size(indexed,1)
                mappings(1,j).Color = Ctmp2(j,:);
                for M = find(indexed(j,:)==1)
                    set(Hlocs(M),'FaceColor',mappings(1,j).Color);
                end
            end
            Nmapping = 1;
            clearvars Ctmp2;
        end
    end
    function import_mapping
        [Mapping_file,Mapping_filepath] = uigetfile('*.xls;*.xlsx','Select File');
        if Mapping_file == 0
            return
        end
        Mapping = lab_read_mappings(fullfile(Mapping_filepath,Mapping_file));
        if isempty(Mapping) | ~isfield(Mapping,'mappingsChannelsFile')
            return
        end
        color = [0.6 0.6 0.6];
        if Mapping.mappingsChannels ~= numchans
            Mapping = lab_reduce_mappings(Mapping);
            if Mapping.mappingsChannels ~= numchans
                return
            end
        end
        indexed = zeros(size(Mapping.mappings,2),numchans);
        for j = 1:size(Hlocs,2)
            set(Hlocs(j),'FaceColor',color);
        end
        mappings = [];
        Nmapping = 1;
        Ctmp2 = get_colormap(size(Mapping.mappings,2));
        for j = 1:size(Mapping.mappings,2)
            if ~isempty(Mapping.mappings{1,j})
                indexed(j,Mapping.mappings{1,j}) = 1;
                mappings(1,j).Name = Mapping.mappingstitle{j,1};
                mappings(1,j).Color = Ctmp2(j,:);
            end
            for M = find(indexed(j,:)==1)
                set(Hlocs(M),'FaceColor',mappings(1,j).Color);
            end
        end
        clearvars Ctmp2;
    end
    function save_mapping
        if isfield(mappings,'active')
            return
        end
        if size(indexed,1) > size(mappings,2)
            indexed = indexed(1:size(mappings,2),:);
        end
        mappingsout.mappingsChannels = numchans;
        mappingsout.mappingsChannelsFile = numchans;
        for j = 1:size(indexed,1)
            mappingsout.mappingstitle{j,1} = mappings(1,j).Name;
            mappingsout.mappings{1,j} = find(indexed(j,:)==1);
        end
        lab_write_mapping([],mappingsout);
    end
    function show_mappinglegend
        if isempty(Hlegend) & max(indexed(:)) > 0
            for j = 1:size(indexed,1)
                Mhandles(j) = find(indexed(j,:)==1,1,'first'); %#ok<AGROW>
                Mnames{j} = mappings(j).Name; %#ok<AGROW>
            end
            Hlegend = legend(Hlocs(Mhandles),Mnames,'Location','EastOutside','FontSize',8);
            Htmp = get(gcf,'Position');
            Htmp(3) = ceil(Htmp(3)*1.3);
            set(gcf,'Position',Htmp);
        elseif ~isempty(Hlegend)
            legend('off');
            Hlegend = [];
            Htmp = get(gcf,'Position');
            Htmp(3) = floor(Htmp(3)/1.3);
            set(gcf,'Position',Htmp);
        end
    end
    function Colormap = get_colormap(Cnumber)
        if ~exist('Cnumber','var')
            if exist('indexed','var') & ~isempty(indexed)
                if size(indexed,1) == 1
                    Cnumber = size(indexed,2);
                else
                    Cnumber = size(indexed,1);
                end
            else
                Colormap = [];
                return
            end
        end
        Colormap = hsv(Cnumber);
        Colormap = Colormap(randperm(Cnumber),:);
    end
    function set_background
        backgroundcolor = get(gcf,'Color');
        if backgroundcolor(1) == 1
            backgroundcolor = [0 0 0];
            textcolor = [1 1 1];
        else
            backgroundcolor = [1 1 1];
            textcolor = [0 0 0];
        end
        set(gcf,'Color',backgroundcolor);
        for T = 1:length(Htext)
            set(Htext(T),'Color',textcolor);
        end
    end
    function edit_scalp
        T = get(Hscalp);
        Prompt = {};
        Formats = {};
        Prompt(end+1,:) = {'Color','FaceColor'};
        Formats(end+1,1).type = 'color';
        Formats(end,1).format = 'color';
        Prompt(end+1,:) = {'Alpha','FaceAlpha'};
        Formats(end+1,1).type = 'range';
        Formats(end,1).limits = [0 1];
        [T,Tskip] = inputsdlg(Prompt,'Set Scalp',Formats,T);
        if Tskip == 0
            set(Hscalp,'FaceAlpha',T.FaceAlpha,'FaceColor',T.FaceColor);
        end
    end
    function edit_brain
        if isempty(Hbrain)
            return
        end
        T = get(Hbrain);
        Prompt = {};
        Formats = {};
        Prompt(end+1,:) = {'Color','FaceColor'};
        Formats(end+1,1).type = 'color';
        Formats(end,1).format = 'color';
        Prompt(end+1,:) = {'Alpha','FaceAlpha'};
        Formats(end+1,1).type = 'range';
        Formats(end,1).limits = [0 1];
        [T,Tskip] = inputsdlg(Prompt,'Set Scalp',Formats,T);
        if Tskip == 0
            set(Hbrain,'FaceAlpha',T.FaceAlpha,'FaceColor',T.FaceColor);
        end
    end
    function view_top
        view(gca,[0 90]); % top
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 0 1])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 0 -1])
            end
        end
    end
    function view_bottom
        view(gca,[180 -90]); % bottom
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 0 -1])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 0 1])
            end
        end
    end
    function view_left
        view(gca,[-90 0]); % left
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[-1 -0.5 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[1 0.5 0])
            end
        end
    end
    function view_right
        view(gca,[90 0]); % right
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[1 -0.5 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[-1 0.5 0])
            end
        end
    end
    function view_front
        view(gca,[180 0]); % front
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 1 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 -1 0])
            end
        end
    end
    function view_back
        view(gca,[0 0]); % back
        if ~isempty(Hlight)
            set(Hlight(1),'Position',[0 -1 0])
            if length(Hlight) > 1
                set(Hlight(2),'Position',[0 1 0])
            end
        end
    end

end

