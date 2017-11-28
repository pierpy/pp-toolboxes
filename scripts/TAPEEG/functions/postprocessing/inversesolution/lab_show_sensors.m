% Plot head-shape and electrodes and store result as picture
%
% lab_show_sensors(LOCS,headshape,filename)
%
% written by F. Hatz 2012

function lab_show_sensors(LOCS,headshape,filename)

if isfield(LOCS,'x')
    elec.chanpos(:,1) = LOCS.x';
    elec.chanpos(:,2) = LOCS.y';
    elec.chanpos(:,3) = LOCS.z';
    elec.elecpos = settings.elec.chanpos;
    elec.label = LOCS.labels';
    elec.unit = 'mm';
    elec = ft_datatype_sens(elec);
elseif isfield(LOCS,'chanpos')
    elec = LOCS;
else
    elec = [];
end

if ~isempty(elec)
    if exist('filename','var')
        f = figure('Visible','off');
    else
        f = figure('Color',[1 1 1],'NumberTitle','off','Name','Sensors','Menubar','none');
        m1 = uimenu(fig1,'Label','File');
        uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
        uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
        uimenu(m1,'Label','Close','Callback','close;');
    end
    axis vis3d; cla
    xlabel('x')
    ylabel('y')
    zlabel('z')
    if exist('headshape','var')
        % plot headshape
        headshape = ft_convert_units(headshape, elec.unit);
        skin = [255 213 119]/255;
        ft_plot_mesh(headshape,'facecolor', skin,'EdgeColor','none','facealpha',0.7);
        lighting gouraud
        material shiny
        camlight
    end
    % plot sensors
    elec = ft_transform_sens(eye(4), elec);
    ft_plot_sens(elec,'label', 'off');
    axis off
    grid on
    view(gca,[-140 15]);
    if exist('filename','var')
        lab_print_figure([filename(1:end-4) '.jpg'],f);
    end
end