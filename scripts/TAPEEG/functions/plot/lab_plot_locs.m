% Plot electrode-locations
%
% [Indexed,Mappings,Locs] = lab_plot_locs(settings,enablesave,editchannels,editmappings,editmontage)
%
%      settings.LOCS                 structure with electrodes-information (see lab_read_locs)
%                                    filepath to electrodesfile
%      settings.Color                color for electrodes
%      settings.indexed              indexed electrodes (integers)
%      settings.ColorIdx             color for indexed electrodes
%      settings.ColorMode            'color' or 'gray'
%      enablesave                    After plotting, script waits for
%                                    save/exit to outpout indexed electrodes
%      editchannels                  1 = editing of channels enabled
%      editmappings                  1 = loading/editing of mappings enabled
%      editmontage                   1 = loading/editing of montages enabled
%
% Written by F. Hatz 2013

function [Indexed,Mappings,Locs] = lab_plot_locs(settings,enablesave,editchannels,editmappings,editmontage,editaux)

if ~exist('editaux','var')
    editaux = 0;
end
if ~exist('editmontage','var')
    if exist('editmappings','var')
        editmontage = editmappings;
    else
        editmontage = 1;
    end
end
if ~exist('editmappings','var')
    editmappings = 1;
end
if ~exist('editchannels','var')
    editchannels = 1;
end
if ~exist('enablesave','var')
    enablesave = 0;
end

Locs = [];
Mappings = [];
Indexed = [];
if ~exist('settings','var') | isempty(settings)
    settings = lab_set_plot_locs;
end
if isfield(settings,'x')
    LOCS = settings;
    clearvars settings
    settings.LOCS = LOCS;
    settings.Color = [1 1 1];
    settings.indexed = [];
    settings.ColorIdx = [1 0 0];
    settings.ColorMode = 'color';
    clearvars LOCS
end
if isempty(settings) | ~isfield(settings,'LOCS') | isempty(settings.LOCS)
    return
end

numchans = size(settings.LOCS.x,2);
NumChanInitial = numchans;

if ~isfield(settings,'Color')
    settings.Color = [1 1 1];
end
if ~isfield(settings,'indexed')
    settings.indexed = [];
end
if ~isfield(settings,'ColorIdx')
    settings.ColorIdx = [0 0 0];
end
if ~isfield(settings,'ColorMode')
    settings.ColorMode = 'color';
end
Color = settings.Color;
if min(Color) == 1
    ColorB = [0 0 0];
else
    ColorB = Color;
end
    
Hlegend = [];

if iscell(settings.indexed) & ~isempty(settings.indexed) & isnumeric(settings.indexed{1})
    Nmapping = 0;
    Ctmp = get_colormap(length(settings.indexed));
    for i = 1:length(settings.indexed)
        Nmapping = Nmapping + 1;
        mappings(1,Nmapping).Name = ['Mapping' num2str(Nmapping)]; %#ok<AGROW>
        mappings(1,Nmapping).Color = Ctmp(Nmapping,:); %#ok<AGROW>
        indexed(Nmapping,:) = zeros(1,numchans); %#ok<AGROW>
        indexed(Nmapping,settings.indexed{Nmapping}) = 1; %#ok<AGROW>
    end
    clearvars Ctmp;
elseif isstruct(settings.indexed) & ~isempty(settings.indexed) & isfield(settings.indexed,'mappings')
    Nmapping = 0;
    Ctmp = get_colormap(length(settings.indexed.mappings));
    for i = 1:length(settings.indexed.mappings)
        Nmapping = Nmapping + 1;
        mappings(1,Nmapping).Name = settings.indexed.mappingstitle{Nmapping,1}; %#ok<AGROW>
        mappings(1,Nmapping).Color = Ctmp(Nmapping,:); %#ok<AGROW>
        indexed(Nmapping,:) = zeros(1,numchans); %#ok<AGROW>
        indexed(Nmapping,settings.indexed.mappings{Nmapping}) = 1; %#ok<AGROW>
    end
    clearvars Ctmp;
else
    Nmapping = 1;
    mappings(1).Name = 'Mapping1';
    mappings(1).Color = settings.ColorIdx;
    indexed = zeros(1,numchans);
    if isnumeric(settings.indexed) & ~isempty(settings.indexed)
        indexed(settings.indexed) = 1;
        Indexed = settings.indexed;
    end
end

settings.PLOT_file = [];
settings.NoMenu = 1;
if isfield(settings,'Title')
    PLOT.Title = settings.Title;
end
PLOT.MinValue = 0;
PLOT.MaxValue = 1;
PLOT.AddPlot = 0;
PLOT.Color = lab_create_cmap(mappings(1).Color,settings.Color);
if isfield(settings,'LOCS_file')
    [~,~,~,PLOT.Title] = lab_filename(settings.LOCS_file);
elseif ~isfield(PLOT,'Title')
    PLOT.Title = 'LOCS';
end

settings = lab_plot_chans(indexed(1,:),PLOT,settings,1);
if ~isfield(settings,'handleL')
    indexed = [];
    return
end
Lhandle = settings.handleL;
Lhandle = Lhandle(:);
Thandle = settings.handleT;
Thandle = Thandle(:);
Chandle = zeros(numchans,numchans);
for i = 1:length(Lhandle)
    set(Thandle(i),'UserData',i,'ButtonDownFcn',@set_channel)
end
if size(indexed,1) > 1
    for m = 1:numchans
        for n = 2:size(indexed,1)
            if indexed(n,m) == 1
                set(Lhandle(m),'FaceColor',mappings(1,n).Color,'EdgeColor',mappings(1,n).Color);
            end
        end
    end
end

if isfield(settings.LOCS,'aux') & settings.LOCS.aux > 0
    numaux = settings.LOCS.aux;
    for Naux = 1:numaux
        draw_dot(Naux)
        set(Thandle(end),'UserData',length(Thandle))
        if editaux == 1
            Ctmp = lines(Naux);
            if isfield(settings.LOCS,'auxNr') & length(settings.LOCS.auxNr) >= Naux & settings.LOCS.auxNr(Naux) > 0
                set(Lhandle(end),'FaceColor',Ctmp(end,:),'EdgeColor',Ctmp(end,:));
                set(Lhandle(settings.LOCS.auxNr(Naux)),'FaceColor',Ctmp(end,:),'EdgeColor',Ctmp(end,:));
            else
                NumAuxInitial = Naux;
            end
        else
            NumAuxInitial = Naux;
        end 
    end
    clearvars Naux
else
    numaux = 0;
    NumAuxInitial = 0;
end

m1 = uimenu(gcf,'Label','File');
uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
if enablesave == 1
    uimenu(m1,'Label','Save','Callback',@(~,~)do_save,'Separator','on');
    uimenu(m1,'Label','Cancel','Callback',@(~,~)do_cancel);
end
m2 = uimenu(gcf,'Label','Load/Save');
uimenu(m2,'Label','Load Indexed','Callback',@(~,~)import_indexed);
uimenu(m2,'Label','Export Indexed','Callback',@(~,~)save_indexed);
if editchannels == 1
    uimenu(m2,'Label','Save LOCS','Callback',@(~,~)save_locs,'Separator','on');
end
m3 = uimenu(gcf,'Label','Edit');
uimenu(m3,'Label','Invert selection','Callback',@(~,~)invert_selection);
if editchannels == 1
    uimenu(m3,'Label','Exclude channels','Callback',@(~,~)del_excluded);
    uimenu(m3,'Label','Exclude selected','Callback',@(~,~)del_selected);
end
uimenu(m3,'Label','Switch ColorMode','Callback',@(~,~)switch_ColorMode,'Separator','on');
uimenu(m3,'Label','Text-Format','Callback',@(~,~)edit_text);
if editmappings == 1
    m4 = uimenu(gcf,'Label','Mapping');
    uimenu(m4,'Label','Select','Callback',@(~,~)select_mapping);
    uimenu(m4,'Label','Edit','Callback',@(~,~)edit_mapping);
    uimenu(m4,'Label','Add','Callback',@(~,~)add_mapping);
    uimenu(m4,'Label','Delete','Callback',@(~,~)del_mapping);
    uimenu(m4,'Label','Legend on/off','Callback',@(~,~)show_mappinglegend,'Separator','on');
    uimenu(m4,'Label','Load','Callback',@(~,~)import_mapping,'Separator','on');
    uimenu(m4,'Label','Save','Callback',@(~,~)save_mapping);
end
if editmontage == 1
    m5 = uimenu(gcf,'Label','Montage');
    uimenu(m5,'Label','Edit','Callback',@(~,~)edit_montage);
    
    uimenu(m5,'Label','Load','Callback',@(~,~)import_montage,'Separator','on');
    uimenu(m5,'Label','Save','Callback',@(~,~)save_montage);
end

if enablesave == 1
    U.Hprint = [];
    U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'normalized',...
        'position', [0.875 0.01 0.12 0.05],...
        'string','Save','fontsize',12,'Callback',@(~,~)do_save,...
        'TooltipString','Save');
    U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'normalized',...
        'position', [0.75 0.01 0.12 0.05],...
        'string','Cancel','fontsize',12,'Callback',@(~,~)do_cancel,...
        'TooltipString','Cancel');
    set(settings.handleF,'UserData',U);
    uiwait
end

    function do_save
        if size(indexed,1) > 1
            for I = 1:size(indexed,1)
                Indexed{I} = find(indexed(I,:)==1);
            end
        else
            Indexed = find(indexed==1);
        end
        Locs = settings.LOCS;
        Mappings = mappings;
        close;
    end
    function do_cancel
        Indexed = [];
        Mappings = [];
        Locs = [];
        close;
    end
    function set_channel(H,~)
        H = get(H);
        T = get(Lhandle(H.UserData));
        if max(T.EdgeColor) == 0
            indexed(Nmapping,H.UserData) = 1;
            if editaux == 0
                set(Lhandle(H.UserData),'FaceColor',mappings(1,Nmapping).Color,'EdgeColor',mappings(1,Nmapping).Color);
            elseif H.UserData > NumChanInitial
                set(Lhandle(H.UserData),'FaceColor',mappings(1,Nmapping).Color,'EdgeColor',mappings(1,Nmapping).Color);
                settings.LOCS.auxNr(H.UserData-NumChanInitial) = H.UserData;
            else
                numaux = numaux + 1;
                tmp = lines(numaux);
                set(Lhandle(H.UserData),'FaceColor',tmp(end,:),'EdgeColor',tmp(end,:));
                auxlabel = inputdlg('test','AUX',[1 20],{['AUX' num2str(numaux)]});
                settings.LOCS.auxlabels{numaux} = auxlabel{1};
                settings.LOCS.auxNr(numaux) = H.UserData;
                settings.LOCS.aux = numaux;
                draw_dot(numaux);
                set(Thandle(end),'UserData',length(Thandle));
                set(Lhandle(end),'FaceColor',tmp(end,:),'EdgeColor',tmp(end,:));
            end
        else
            indexed(:,:) = 0;
            set(Lhandle,'EdgeColor',ColorB,'FaceColor',Color);
            if editaux == 1
                if length(Thandle) > NumAuxInitial+NumChanInitial
                    delete(Thandle(NumAuxInitial+NumChanInitial+1:end));
                    delete(Lhandle(NumAuxInitial+NumChanInitial+1:end));
                    Thandle = Thandle(1:NumAuxInitial+NumChanInitial);
                    Lhandle = Lhandle(1:NumAuxInitial+NumChanInitial);
                    Chandle = Chandle(1:NumAuxInitial+NumChanInitial,1:NumAuxInitial+NumChanInitial);
                    numchans = length(Thandle);
                    if isfield(settings.LOCS,'auxlabels') & ~isempty(settings.LOCS.auxlabels)
                        if NumAuxInitial == 0
                            settings.LOCS.auxlabels = {};
                        else
                            settings.LOCS.auxlabels = settings.LOCS.auxlabels(1:NumAuxInitial);
                        end
                    end
                    settings.LOCS.auxNr = [];
                    numaux = NumAuxInitial;
                    settings.LOCS.aux = NumAuxInitial;
                end
            end
        end
    end
    function import_indexed
        if isfield(mappings,'active')
            indexed = zeros(1,numchans);
            for j = 1:length(Lhandle)
                set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color);
            end
            delete(Chandle(Chandle>0));
            Chandle = zeros(numchans,numchans);
            mappings = [];
            mappings(1).Name = 'Mapping1';
            mappings(1).Color = settings.ColorIdx;
            Nmapping = 1;
        end
        [INDEX_file,INDEX_filepath] = uigetfile('*.xls;*.xlsx','Select File');
        if INDEX_file == 0
            return
        end
        if ispc
            [nchans,~,xlsinput] = xlsread(fullfile(INDEX_filepath,INDEX_file),1);
        else
            [nchans,~,xlsinput] = xlsread(fullfile(INDEX_filepath,INDEX_file),1,'','basic');
        end
        if ~isempty(find(nchans==numchans-numaux,1)) & size(xlsinput,2) == 2
            nchans = xlsinput{find(nchans==numchans-numaux,1),2};
            if ischar(nchans)
                nchans = str2num(nchans); %#ok<ST2NM>
            end
        elseif size(xlsinput,2) == 2 & ischar(xlsinput{1,2})
            nchans = [];
        elseif size(nchans,1) > 1 & size(nchans,2) == 1
            nchans = nchans';
        elseif size(nchans,1) > 1
            nchans = nchans(1,:);
        end
        for j = 1:numchans-numaux
            if indexed(Nmapping,j) == 1
                set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color);
            end
        end
        indexed(Nmapping,:) = zeros(1,numchans-numaux);
        if ~isempty(nchans)
            indexed(Nmapping,nchans) = 1;
        end
        for j = 1:numchans-numaux
            if indexed(Nmapping,j) == 1
                set(Lhandle(j),'FaceColor',mappings(1,Nmapping).Color,'EdgeColor',mappings(1,Nmapping).Color);
            end
        end
    end
    function save_indexed
        [INDEX_file,INDEX_filepath] = uiputfile('*.xls;*.xlsx','Select File to store');
        if INDEX_file == 0
            return
        end
        if size(mappings,2) > 1
            [~,~,format,INDEX_file] = lab_filename(INDEX_file);
            INDEX_file = [INDEX_file '_' mappings(1,Nmapping).Name '.' format];
            clearvars format
        end
        xlsout{1,1} = numchans;
        xlsout{1,2} = num2str(find(indexed(Nmapping,:)==1));
        lab_write_xls(fullfile(INDEX_filepath,INDEX_file),xlsout);
    end
    function save_locs
        [LOCS_file,LOCS_filepath] = uiputfile('*.els','Select File to store');
        if LOCS_file == 0
            return
        end
        if exist(fullfile(LOCS_filepath,LOCS_file),'file')
            delete(fullfile(LOCS_filepath,LOCS_file));
        end
        lab_write_els(fullfile(LOCS_filepath,LOCS_file),settings.LOCS);
    end
    function del_excluded
        header = [];
        if isfield (settings,'LOCS') & isfield(settings.LOCS,'labels')
            header.channels = char(settings.LOCS.labels);
            header.numchannels = size(header.channels,1);
        end
        settings = lab_set_exclude(settings,header);
        if isfield(settings,'exclude') & ~isempty(settings.exclude) 
            for j = 1:length(settings.exclude)
                delete(Lhandle(settings.exclude(j)));
                delete(Thandle(settings.exclude(j)));
            end
            includechannels = 1:numchans;
            includechannels = setdiff(includechannels,settings.exclude);
            settings.LOCS.x = settings.LOCS.x(1,includechannels);
            settings.LOCS.y = settings.LOCS.y(1,includechannels);
            settings.LOCS.z = settings.LOCS.z(1,includechannels);
            settings.LOCS.labels = settings.LOCS.labels(1,includechannels);
            if isfield(settings.LOCS,'radius')
                settings.LOCS.radius = settings.LOCS.radius(1,includechannels);
                settings.LOCS.theta = settings.LOCS.theta(1,includechannels);
            end
            if isfield(settings.LOCS,'sph_radius')
                settings.LOCS.sph_radius = settings.LOCS.sph_radius(1,includechannels);
                settings.LOCS.sph_theta = settings.LOCS.sph_theta(1,includechannels);
                settings.LOCS.sph_phi = settings.LOCS.sph_phi(1,includechannels);
            end
            indexed = indexed(:,includechannels);
            Thandle = Thandle(includechannels);
            Lhandle = Lhandle(includechannels);
            numchans = length(includechannels);
            for j = 1:length(Lhandle)
                set(Thandle(j),'UserData',j)
            end
            settings = rmfield(settings,'exclude');
        end
    end
    function del_selected
        selected = find(indexed(Nmapping,:) == 1);
        if ~isempty(selected)
            for j = 1:length(selected)
                delete(Lhandle(selected(j)));
                delete(Thandle(selected(j)));
            end
            includechannels = 1:numchans;
            includechannels = setdiff(includechannels,selected);
            settings.LOCS.x = settings.LOCS.x(1,includechannels);
            settings.LOCS.y = settings.LOCS.y(1,includechannels);
            settings.LOCS.z = settings.LOCS.z(1,includechannels);
            settings.LOCS.labels = settings.LOCS.labels(1,includechannels);
            if isfield(settings.LOCS,'radius')
                settings.LOCS.radius = settings.LOCS.radius(1,includechannels);
                settings.LOCS.theta = settings.LOCS.theta(1,includechannels);
            end
            indexed = indexed(:,includechannels);
            Thandle = Thandle(includechannels);
            Lhandle = Lhandle(includechannels);
            numchans = length(includechannels);
            for j = 1:length(Lhandle)
                set(Thandle(j),'UserData',j)
            end
            if size(mappings,2) == 1
                mappings.Name = 'Mapping1';
                indexed = zeros(1,numchans);
                for j = 1:length(Lhandle)
                    set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color);
                end
            else
                include = setdiff(1:size(mappings,2),Nmapping);
                mappings = mappings(1,include);
                indexed = indexed(include,:);
                Nmapping = 1;
            end
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
                set(Lhandle(M),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
            end
        end
        clearvars Ctmp2;
        mappings(1,Nmapping) = inputsdlg(mappings(1,Nmapping),'Mappings');
    end
    function edit_mapping
        if isfield(mappings,'active')
            return
        end
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
                set(Lhandle(M),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
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
            for j = 1:length(Lhandle)
                set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color,'LineWidth',0.5);
            end
        else
            include = setdiff(1:size(mappings,2),Nmapping);
            for j = find(indexed(Nmapping,:)==1)
                set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color);
            end
            mappings = mappings(1,include);
            indexed = indexed(include,:);
            Ctmp2 = get_colormap(size(mappings,2));
            for j = 1:size(indexed,1)
                mappings(1,j).Color = Ctmp2(j,:);
                for M = find(indexed(j,:)==1)
                    set(Lhandle(M),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
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
        if Mapping.mappingsChannels ~= numchans-numaux
            Mapping = lab_reduce_mappings(Mapping);
            if Mapping.mappingsChannels ~= numchans-numaux
                return
            end
        end
        indexed = zeros(size(Mapping.mappings,2),numchans);
        for j = 1:length(Lhandle)
            set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color,'LineWidth',0.5);
        end
        delete(Chandle(Chandle>0));
        Chandle = zeros(numchans,numchans);
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
                set(Lhandle(M),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
            end
        end
        clearvars Ctmp2;
    end
    function save_mapping
        if isfield(mappings,'active')
            return
        end
        mappingsout.mappingsChannels = numchans-numaux;
        mappingsout.mappingsChannelsFile = numchans-numaux;
        for j = 1:size(indexed,1) - numaux
            mappingsout.mappingstitle{j,1} = mappings(1,j).Name;
            mappingsout.mappings{1,j} = find(indexed(j,:)==1);
        end
        lab_write_mapping([],mappingsout);
    end
    function show_mappinglegend
        if isempty(Hlegend) & max(indexed(:)) > 0 & ~isfield(mappings,'active')
            for j = 1:size(indexed,1)
                Mhandles(j) = find(indexed(j,:)==1,1,'first'); %#ok<AGROW>
                Mnames{j} = mappings(j).Name; %#ok<AGROW>
            end
            Hlegend = legend(Lhandle(Mhandles),Mnames,'Location','EastOutside','FontSize',8);
            tmp = get(gcf,'Position');
            tmp(3) = ceil(tmp(3)*1.3);
            set(gcf,'Position',tmp);
        elseif ~isempty(Hlegend)
            legend('off');
            Hlegend = [];
            tmp = get(gcf,'Position');
            tmp(3) = floor(tmp(3)/1.3);
            set(gcf,'Position',tmp);
        end
    end
    function invert_selection
        if size(indexed,1) > 1
            return
        end
        Tindexed = zeros(1,numchans);
        Tindexed(indexed==0) = 1;
        indexed = Tindexed;
        clearvars Tindexed;
        for j = 1:numchans
            if indexed(1,j) == 1
                set(Lhandle(j),'FaceColor',mappings(1,1).Color,'EdgeColor',mappings(1,1).Color,'LineWidth',0.5);
            else
                set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color,'LineWidth',0.5);
            end
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
        if isfield(settings,'ColorMode') & strcmp(settings.ColorMode,'color')
            Colormap = [];
            nmap = Cnumber + 1;
            while size(Colormap,1) < Cnumber
                Colormap = AdvancedColormap('wkrwmgwcbwy',nmap);
                Colormap = Colormap(sum(Colormap,2)<2.9,:);
                Colormap = Colormap(sum(Colormap,2)>0.45,:);
                nmap = nmap + 1;
            end
            clearvars nmap
            Colormap = Colormap(randperm(Cnumber),:);
        else
            Colormap = flipud(gray(round((Cnumber+1)*1.3333)));
            Colormap = Colormap(2:Cnumber+1,:);
            Colormap = Colormap(randperm(Cnumber),:);
        end
    end

    function switch_ColorMode(H1,H2)
        if isfield(settings,'ColorMode') & strcmp(settings.ColorMode,'color')
            settings.ColorMode = 'gray';
        else
            settings.ColorMode = 'color';
        end
        Color = [1 1 1];
        ColorB = [0 0 0];
        Ctmp2 = get_colormap(size(mappings,2));
        for j = 1:size(indexed,1)
            if ~isfield(mappings,'active')
                mappings(1,j).Color = Ctmp2(j,:);
            else
                if strcmp(settings.ColorMode,'gray')
                    mappings(1,j).Color = [0.7 0.7 0.7];
                else
                    mappings(1,j).Color = [0 1 0];
                end
            end
            for M = find(indexed(j,:)==1)
                if isempty(find(indexed(:,M) == -1,1))
                    set(Lhandle(M),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
                else
                    set(Lhandle(M),'FaceColor',mappings(1,j).Color);
                end
            end
            if isfield(mappings,'active')
                for M = find(indexed(j,:) == -1)
                    if strcmp(settings.ColorMode,'gray')
                        set(Lhandle(M),'EdgeColor',[0.4 0.4 0.4],'LineWidth',2);
                    else
                        set(Lhandle(M),'EdgeColor',[1 0 0],'LineWidth',2);
                    end
                end
                for M = find(indexed(j,:) == 1)
                    if strcmp(mappings(1,j).reference,'LAPL')
                        if strcmp(settings.ColorMode,'gray')
                            set(Lhandle(M),'EdgeColor',[0.2 0.2 0.2],'LineWidth',2);
                        else
                            set(Lhandle(M),'EdgeColor',[0 0 1],'LineWidth',2);
                        end
                    end
                end
            end
        end
        tmp = find(sum(indexed,1)==0);
        for j = 1:length(tmp)
            set(Lhandle(tmp(j)),'EdgeColor',ColorB,'FaceColor',Color,'LineWidth',0.5);
        end
        tmp = find(Chandle>0);
        if strcmp(settings.ColorMode,'gray');
            for j = tmp
                set(Chandle(j),'Color',[0.3 0.3 0.3]);
            end
        else
            for j = tmp
                set(Chandle(j),'Color',[1 0 0]);
            end
        end
    end

    function import_montage
        [Montage_file,Montage_filepath] = uigetfile('*.xls;*.xlsx','Select File');
        if Montage_file == 0
            return
        end
        Montage = lab_read_montage(fullfile(Montage_filepath,Montage_file));
        if isempty(Montage) | ~isfield(Montage,'chans')
            return
        elseif ~isempty(Montage(1,1).numchans) & Montage(1,1).numchans ~= numchans-numaux
            cfg.EXTRA.numdatachans = numchans-numaux;
            Montage = lab_reduce_montage(Montage,cfg);
            if Montage(1,1).numchans ~= numchans-numaux
                disp('    Abort importing Montage, mismatch in number of channels')
                return
            end
        end
        if size(Montage,2) > 1
            for j = 1:size(Montage,2)
                strlist{j,1} = Montage(1,j).name; %#ok<AGROW>
            end
            selection = listdlg('PromptString','Select Montage','SelectionMode','single','ListString',strlist);
            if isempty(selection)
                return
            else
                Montage = Montage(1,selection);
            end
            clearvars strlist selection
        end
        plot_montage(Montage);
    end
    function save_montage
        if ~isfield(mappings,'active') | ~isfield(mappings,'reference')
            return
        end
        MontageOut.numchans = numchans-numaux;
        if isfield(mappings,'Filename')
            MontageOut.name = mappings(1,1).Filename;
        else
            MontageOut.name = 'Montage';
        end
        for j = 1:size(indexed,1)
            MontageOut.label{j,1} = mappings(1,j).Name;
            MontageOut.chans{j,1} = find(indexed(j,:)==1);
            MontageOut.chans{j,2} = 0;
            if max(MontageOut.chans{j,1}) > (numchans-numaux)
                MontageOut.chans{j,1} = MontageOut.chans{j,1}(1) - (numchans-numaux);
                MontageOut.chans{j,2} = 1;
            end
            if length(mappings(1,j).reference) >= 3 & strcmp(mappings(1,j).reference(1:3),'AVG')
                MontageOut.chans{j,3} = 'AVG';
                MontageOut.chans{j,4} = 0;
            elseif length(mappings(1,j).reference) > 3 & strcmp(mappings(1,j).reference(1:4),'LAPL')
                MontageOut.chans{j,3} = 'LAPL';
                MontageOut.chans{j,4} = 0;
            else
                MontageOut.chans{j,3} = find(indexed(j,:)==-1);
                MontageOut.chans{j,4} = 0;
                if max(MontageOut.chans{j,3}) > (numchans-numaux)
                    MontageOut.chans{j,3} = MontageOut.chans{j,3}(1) - (numchans-numaux);
                    MontageOut.chans{j,4} = 1;
                end
            end
        end
        lab_write_montage([],MontageOut);
    end

    function edit_montage
        if ~isfield(mappings,'active') | ~isfield(mappings,'reference')
            return
        end
        if isfield(settings.LOCS,'auxlabels')
            auxlabels = settings.LOCS.auxlabels;
        else
            auxlabels = {};
        end
        for j = length(auxlabels)+1:6
            auxlabels{1,j} = ['AUX' num2str(j)];
        end
        ListActive = [settings.LOCS.labels auxlabels];
        ListRef = [settings.LOCS.labels auxlabels {'AVG','LAPL'}];
        for j = 1:size(indexed,1)
            MontageList{j,1} = mappings(1,j).Name; %#ok<AGROW>
            
            if length(find(indexed(j,:)==1)) == 1
                MontageList{j,2} = ListActive{indexed(j,:)==1}; %#ok<AGROW>
            elseif  length(find(indexed(j,:)==1)) > 1
                tmp = find(indexed(j,:)==1);
                MontageList{j,2} = ''; %#ok<AGROW>
                for L = tmp
                    MontageList{j,2} = [MontageList{j,2} ListActive{L} ' ']; %#ok<AGROW>
                end
                if max(strcmp(ListRef,MontageList{j,2})) == 0
                    ListActive = [ListActive MontageList(j,2)]; %#ok<AGROW>
                end
            end
            
            if length(mappings(1,j).reference) >= 3 & strcmp(mappings(1,j).reference(1:3),'AVG')
                MontageList{j,3} = mappings(1,j).reference; %#ok<AGROW>
            elseif length(mappings(1,j).reference) > 3 & strcmp(mappings(1,j).reference(1:4),'LAPL')
                MontageList{j,3} = mappings(1,j).reference; %#ok<AGROW>
            elseif length(find(indexed(j,:)==-1)) == 1
                MontageList{j,3} = ListRef{indexed(j,:)==-1}; %#ok<AGROW>
            elseif  length(find(indexed(j,:)==-1)) > 1
                tmp = find(indexed(j,:)==-1);
                MontageList{j,3} = ''; %#ok<AGROW>
                for L = tmp
                    MontageList{j,3} = [MontageList{j,3} ListRef{L} ' ']; %#ok<AGROW>
                end
                if max(strcmp(ListRef,MontageList{j,3})) == 0
                    ListRef = [ListRef MontageList(j,3)]; %#ok<AGROW>
                end
            end
        end
        MontageList = lab_table_dialog(MontageList,{'Name','active','reference'},'Edit Montage',2, ...
            {'char',ListActive,ListRef});
        if isempty(MontageList)
            return
        else
            for j = 1:size(MontageList,1)
                for M = 2:3
                    tmp1 = textscan(MontageList{j,M},'%s ');
                    tmp1 = tmp1{1,1};
                    tmp2 = [];
                    for N = 1:size(tmp1,1)
                        if ~isempty(find(strcmp(ListActive,tmp1{N,1}),1))
                            tmp2 = [tmp2 find(strcmp(ListActive,tmp1{N,1}),1)]; %#ok<AGROW>
                        end
                    end
                    if ~isempty(tmp2)
                        MontageList{j,M} = tmp2;
                    end
                end
            end
            clearvars j M N tmp1 tmp2
        end
        for j = 1:size(MontageList,1)
            if ischar(MontageList{j,2}) | max(MontageList{j,2}) <= numchans - numaux
                Montage.chans{j,1} = MontageList{j,2};
                Montage.chans{j,2} = 0;
            elseif isnumeric(MontageList{j,2})
                tmp = MontageList{j,2};
                tmp = tmp - (numchans-numaux);
                tmp = tmp(tmp>0);
                Montage.chans{j,1} = tmp;
                Montage.chans{j,2} = 1;
            end
            if ischar(MontageList{j,3}) | max(MontageList{j,3}) <= numchans - numaux
                Montage.chans{j,3} = MontageList{j,3};
                Montage.chans{j,4} = 0;
            elseif isnumeric(MontageList{j,3})
                tmp = MontageList{j,3};
                tmp = tmp - (numchans-numaux);
                tmp = tmp(tmp>0);
                Montage.chans{j,3} = tmp;
                Montage.chans{j,4} = 1;
            end
        end
        Montage.label = MontageList(:,1);
        if isfield(mappings,'Filename')
            Montage.name = mappings(1,1).Filename;
        else
            Montage.name = 'Montage';
        end
        plot_montage(Montage);
    end
    
    function plot_montage(Montage)
        indexed = zeros(size(Montage.chans,2),numchans);
        for j = 1:length(Lhandle)
            set(Lhandle(j),'EdgeColor',ColorB,'FaceColor',Color,'LineWidth',0.5);
        end
        delete(Chandle(Chandle>0));
        Chandle = zeros(numchans,numchans);
        mappings = [];
        Nmapping = 1;
        mappings.Filename = Montage.name;
        if strcmp(settings.ColorMode,'gray')
            edgecolor = [0.4 0.4 0.4];
            edgecolor2 = [0.2 0.2 0.2];
        else
            edgecolor = [1 0 0];
            edgecolor2 = [0 0 1];
        end
        for j = 1:size(Montage.chans,1)
            if strcmp(settings.ColorMode,'gray')
                mappings(1,j).Color = [0.7 0.7 0.7];
            else
                mappings(1,j).Color = [0 1 0];
            end
            if Montage.chans{j,2} == 0 & min(Montage.chans{j,1}) >= 1 & max(Montage.chans{j,1}) <= numchans-numaux
                indexed(j,Montage.chans{j,1}) = 1;
                mappings(1,j).active = '';
                for M = find(indexed(j,:)==1)
                    if isempty(find(indexed(:,M) == -1,1))
                        set(Lhandle(M),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
                    else
                        set(Lhandle(M),'FaceColor',mappings(1,j).Color);
                    end
                end
            elseif Montage.chans{j,2} == 1 & length(Montage.chans{j,1}) == 1 & Montage.chans{j,1} >= 1 & Montage.chans{j,1} <= 6
                if numaux < Montage.chans{j,1}
                    for N = numaux+1:Montage.chans{j,1}
                        draw_dot(N);
                    end
                    numaux = Montage.chans{j,1};
                end
                indexed(j,numchans-numaux+Montage.chans{j,1}) = 1;
                set(Lhandle(indexed(j,:)==1),'FaceColor',mappings(1,j).Color,'EdgeColor',mappings(1,j).Color);
                mappings(1,j).active = '';
            end
            
            if Montage.chans{j,4} == 0 & isnumeric(Montage.chans{j,3}) & ...
                    min(Montage.chans{j,3}) >= 1 & max(Montage.chans{j,3}) <= numchans
                tmp = Montage.chans{j,3};
                if Montage.chans{j,2} == 0 & min(Montage.chans{j,1}) >= 1
                    tmp = setdiff(tmp,Montage.chans{j,1});
                end
                indexed(j,tmp) = -1;
                mappings(1,j).reference = '';
                for M = find(indexed(j,:) == -1)
                    set(Lhandle(M),'EdgeColor',edgecolor,'LineWidth',2);
                end
            elseif Montage.chans{j,4} == 1 & length(Montage.chans{j,3}) == 1 & Montage.chans{j,3} >= 1 & Montage.chans{j,3} <= 6
                if numaux < Montage.chans{j,3}
                    for N = numaux+1:Montage.chans{j,3}
                        draw_dot(N);
                    end
                    numaux = Montage.chans{j,3};
                end
                indexed(j,numchans-numaux+Montage.chans{j,3}) = -1;
                mappings(1,j).reference = '';
                for M = find(indexed(j,:) == -1)
                    set(Lhandle(M),'EdgeColor',edgecolor,'LineWidth',2);
                end
            elseif Montage.chans{j,4} == 0 & ischar(Montage.chans{j,3})
                mappings(1,j).reference = Montage.chans{j,3};
                if strcmp(Montage.chans{j,3},'LAPL')
                    for M = find(indexed(j,:) == 1)
                        set(Lhandle(M),'EdgeColor',edgecolor2,'LineWidth',2);
                    end
                elseif strcmp(Montage.chans{j,3},'AVG')
                    for M = find(indexed(j,:) == 1)
                        set(Lhandle(M),'EdgeColor',[0 0 0],'LineWidth',2);
                    end
                end
            end
            mappings(1,j).Name = Montage.label{j,1};
            if ~isempty(find(indexed(j,:)==1,1)) & ~isempty(find(indexed(j,:)==-1,1))
                P1 = find(indexed(j,:)==1);
                P2 = find(indexed(j,:)==-1);
                for p1 = P1
                    for p2 = P2
                        if strcmp(settings.ColorMode,'gray')
                            lab_drawline(p1,p2,[0.3 0.3 0.3],1.5);
                        else
                            lab_drawline(p1,p2,[1 0 0],1.5);
                        end
                    end
                end
            end
        end
        clearvars Ctmp2;
    end

    function lab_drawline(P1,P2,Color,Size)
        if ~exist('Size','var')
            Size = 1;
        end
        if ~exist('Color','var')
            Color = [0 0 0];
        end
        X(1,1) = mean(get(Lhandle(P1),'XData'));
        X(1,2) = mean(get(Lhandle(P2),'XData'));
        Y(1,1) = mean(get(Lhandle(P1),'YData'));
        Y(1,2) = mean(get(Lhandle(P2),'YData'));
        Chandle(P1,P2) = line(X,Y,[-0.0005 -0.0005],'LineWidth',Size,'Color',Color);
    end

    function draw_dot(N)
        Lhandle(end+1,1) = copyobj(Lhandle(end),gca);
        Thandle(end+1,1) = copyobj(Thandle(end),gca);
        Chandle(end+1,end+1) = 0;
        numchans = numchans+1;
        Lpos = get(Lhandle(end),'Vertices');
        Pos(1) = mean(Lpos(:,1));
        Pos(2) = mean(Lpos(:,2));
        Lpos(:,1) = Lpos(:,1) - Pos(1) - 0.45 + (N-1)*0.05;
        Lpos(:,2) = Lpos(:,2) - Pos(2) - 0.45;
        set(Lhandle(end),'Vertices',Lpos,'EdgeColor',ColorB,'FaceColor',Color)
        Tpos = get(Thandle(end),'Position');
        Tpos(1) = - 0.45 + (N-1)*0.05;
        Tpos(2) = - 0.45;
        if isfield(settings.LOCS,'auxlabels') & ~isempty(settings.LOCS.auxlabels) & ~isempty(settings.LOCS.auxlabels{N})
            set(Thandle(end),'Position',Tpos,'String',settings.LOCS.auxlabels{N},'HorizontalAlignment','center','VerticalAlignment','middle');
        else
            set(Thandle(end),'Position',Tpos,'String',['AUX' num2str(N)],'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    
    function edit_text
        if ~isempty(Thandle)
            tmp = inputdlg({'Text Size(%)','Text Position'},'Text',[1 15],{'100','0'});
            if ~isempty(tmp)
                fsize = str2double(tmp{1});
                fposition = str2double(tmp{2});
            else
                return
            end
            if ~isempty(fsize)
                for t = 1:length(Thandle)
                    pos = get(Thandle(t),'Position');
                    siz = get(Thandle(t),'FontSize');
                    pos(2) = pos(2) + fposition*siz;
                    siz = siz * fsize/100;
                    set(Thandle(t),'FontSize',siz,'Position',pos);
                end
            end
        end
    end
end