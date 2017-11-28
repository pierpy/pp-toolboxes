function settings = lab_plot_legend(Legend,settings,doedit)

if ~exist('doedit','var')
    doedit = false;
end
if ~exist('settings','var')
    settings = [];
end
if isfield(settings,'Fhandle') & ~isempty(settings.Fhandle) & ishandle(settings.Fhandle)
    figure(settings.Fhandle);
else
    settings.Fhandle = gcf;
end
if ~isfield(settings,'Ahandle') | isempty(settings.Ahandle) | ~ishandle(settings.Ahandle)
    settings.Ahandle = gca;
end
if isfield(settings,'Lhandle') & ~isempty(settings.Lhandle) & ishandle(settings.Lhandle)
    U = get(settings.Lhandle,'UserData');
else
    U.Lhandle = [];
    U.Color = {};
    U.Text = {};
    U.Mode = '';
end

if ischar(Legend) & strcmp(Legend,'off')
    colorbar off
    legend off
    if isfield(U,'PosL') & isfield(U,'PosE')
        set(settings.Fhandle,'Units','pixels');
        set(settings.Fhandle,'Position',U.PosL);
        set(settings.Ahandle,'Units','pixels');
        set(settings.Fhandle,'Position',U.PosE);
        set(settings.Fhandle,'Units','normalized');
        set(settings.Ahandle,'Units','normalized');
    end
    settings.Lhandle = [];
    return
end

hold on
for i = 1:length(Legend)
    if isfield(Legend,'Handle') & ishandle(Legend(i).Handle)
        if isfield(Legend,'Text') & ~isempty(Legend(i).Text)
            if ischar(Legend(i).Text)
                Text = {Legend(i).Text};
            else
                Text = Legend(i).Text;
            end
        else
            Text = {['Object ' num2str(length(U.Lhandle)+1)]};
        end
        [U,Flag] = combine_legend(U,Legend(i).Color,Text);
        if Flag == true
            U.Lhandle(end+1,1) = Legend(i).Handle;
        end
    elseif isfield(Legend,'Color') & ~isempty(Legend(i).Color)
        if isfield(Legend,'Mode') & strcmp(Legend(i).Mode,'Colorbar')
            if isfield(Legend,'Text') & ~isempty(Legend(i).Text)
                Text = Legend(i).Text;
            else
                Text = cellstr(num2str(caxis'));
            end
            U = combine_legend(U,Legend(i).Color,Text);
            clearvars Text
        elseif isfield(Legend,'Mode') & strcmp(Legend(i).Mode,'Line')
            if isfield(Legend,'Text') & ~isempty(Legend(i).Text)
                if ischar(Legend(i).Text)
                    Text = {Legend(i).Text};
                else
                    Text = Legend(i).Text;
                end
            else
                Text = {['Object ' num2str(length(U.Lhandle)+1)]};
            end
            [U,Flag] = combine_legend(U,Legend(i).Color,Text);
            if Flag == 1
                U.Lhandle(end+1,1) = line([0;1.0e-12],[0;1.0e-12],'Color',Legend(i).Color,'Visible','off');
            end
        else
            if isfield(Legend,'Text') & ~isempty(Legend(i).Text)
                if ischar(Legend(i).Text)
                    Text = {Legend(i).Text};
                else
                    Text = Legend(i).Text;
                end
            else
                Text = {['Object ' num2str(length(U.Lhandle)+1)]};
            end
            [U,Flag] = combine_legend(U,Legend(i).Color,Text);
            if Flag == 1
                U.Lhandle(end+1,1) = fill([0 0 1.0e-12],[0 1.0e-12 1.0e-12],Legend(i).Color,'Visible','off');
            end
        end
    end
end
if ~isempty(U.Lhandle)
    set(settings.Fhandle,'Units','pixels');
    set(settings.Ahandle,'Units','pixels');
    if ~isfield(settings,'Lhandle') | isempty(settings.Lhandle) | ~ishandle(settings.Lhandle)
        doresize = true;
    else
        doresize = false;
        legend off
    end
    Text = {};
    for i = 1:length(U.Text)
        Text = cat(1,Text,U.Text{i});
    end
    if length(U.Lhandle) > 12
        FontSize = 7;
    else
        FontSize = 9;
    end
    if doedit == true
        Text = lab_table_dialog(Text,{'Name'},'Names',0);
    end
    settings.Lhandle = legend(U.Lhandle,Text,'Fontsize',FontSize,'Location','East');
    U.Mode = 'Legend';
    set(settings.Lhandle,'Units','pixels');
    Fpos = get(settings.Fhandle,'Position');
    Lpos = get(settings.Lhandle,'Position');
    if doresize == true
        U.PosE = Fpos;
        Lpos(1) = Fpos(3) + 10;
        Lpos(2) = round((Fpos(4)-Lpos(4))/2);
        Fpos(3) = Fpos(3) + Lpos(3) + 20;
        U.PosL = Fpos;
        set(settings.Fhandle,'Position',Fpos);
    else
        Lpos(1) = Fpos(3) - Lpos(3) -10;
        Lpos(2) = round((Fpos(4)-Lpos(4))/2);
    end
    set(settings.Lhandle,'Position',Lpos);
    set(settings.Fhandle,'Units','normalized');
    set(settings.Ahandle,'Units','normalized');
    set(settings.Lhandle,'Units','normalized','UserData',U);
elseif ~isempty(U.Color)
    set(settings.Fhandle,'Units','pixels');
    set(settings.Ahandle,'Units','pixels');
    if ~isfield(settings,'Lhandle') | isempty(settings.Lhandle) | ~ishandle(settings.Lhandle)
        doresize = true;
    else
        doresize = false;
        colorbar off
    end
    Color = [];
    Text = {};
    for i = 1:length(U.Color)
        Color = cat(1,Color,U.Color{i});
        if isempty(Text)
            Text = U.Text{i};
        elseif length(U.Text{i}) > 1
            Text{end} = [U.Text{i}{1} ' | ' Text{end}];
            Text = cat(1,Text,U.Text{i}(2:end));
        end
    end
    set(gcf,'Colormap',Color);
    settings.Lhandle = colorbar('Location','East');
    tmp = get(settings.Lhandle,'YLim');
    tmp = tmp(1):(tmp(2)-tmp(1))/(length(Text)-1):tmp(2);
    set(settings.Lhandle,'YTickLabel',Text,'YTick',tmp(:),'Fontsize',8);
    U.Mode = 'Colorbar';
    set(settings.Lhandle,'Units','pixels');
    Fpos = get(settings.Fhandle,'Position');
    Lpos = get(settings.Lhandle,'Position');
    if doresize == true
        U.PosE = Fpos;
        Lpos(1) = Fpos(3) + 10;
        Lpos(2) = round((Fpos(4)-Lpos(4))/2);
        Fpos(3) = Fpos(3) + Lpos(3) + 20;
        U.PosL = Fpos;
        set(settings.Fhandle,'Position',Fpos);
    else
        Lpos(1) = Fpos(3) - Lpos(3) -10;
        Lpos(2) = round((Fpos(4)-Lpos(4))/2);
    end
    set(settings.Lhandle,'Position',Lpos);
    set(settings.Fhandle,'Units','normalized');
    set(settings.Ahandle,'Units','normalized');
    set(settings.Lhandle,'Units','normalized','UserData',U);
end

end

function [U,Flag] = combine_legend(U,Color,Text)
    Flag = true;
    for j = 1:length(U.Color)
        if max(abs(U.Color{j}(:)-Color(:))) == 0
            for m = 1:length(U.Text{j})
                if m <= length(Text) & strcmp(U.Text{j}{m},Text{m})
                    Flag = false;
                end
            end
        end
    end
    if Flag == true
        U.Color{end+1,1} = Color;
        U.Text{end+1,1} = Text(:);
    end
end