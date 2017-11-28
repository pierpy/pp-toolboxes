% Helper function for lab_collect_spectras
%
% result = lab_collect_spectraldata_plot(cfg,calc,result,dowriting)
%
% written by F. Hatz 2012

function result = lab_collect_spectraldata_plot(cfg,calc,result,dowriting)

if ~exist('dowriting','var')
    dowriting = true;
end

if ~isfield(result,'T') | isempty(result.T)
    return
end

if isfield(cfg.CollectFFT,'correctpf') & cfg.CollectFFT.correctpf == false & dowriting == true
    skipfig = 1;
else
    skipfig = 0;
end
if isfield(cfg.CollectFFT,'correctpf') & cfg.CollectFFT.correctpf == true
    skipedit = 0;
else
    skipedit = 1;
end

f = figure('Visible','off','Color',[1 1 1],'MenuBar','none','NumberTitle','off','Name','Spectrum');

Nplot = 1;
MaxPlot = size(result.T,2);
flagexit = 0;
while flagexit == 0
    m1 = uimenu(f,'Label','File');
    uimenu(m1,'Label','Save Picture','Callback',@(~,~)lab_print_figure);
    uimenu(m1,'Label','Copy Picture','Callback','print -dmeta;');
    uimenu(m1,'Label','Close','Callback',@go_end);
    m2 = uimenu(f,'Label','Edit');
    uimenu(m2,'Label','Previous','Callback',@go_previous);
    uimenu(m2,'Label','Next','Callback',@go_next);
    uimenu(m2,'Label','GoTo','Callback',@go_select);
    if skipedit == 0
        m3 = uimenu(f,'Label','Extra');
        uimenu(m3,'Label','Set NaN','Callback',@(hObj,evd)set_NaN(hObj,evd,cfg.CollectFFT));
        uimenu(m3,'Label','Set Limits','Callback',@(hObj,evd)set_Limits(hObj,evd),'Separator','on');
    end
    PLOT.T = result.T(1,Nplot);
    PLOT.R = result.R(1,Nplot);
    PLOT.patient = result.patient{1,Nplot};
    if cfg.CollectFFT.qualityplot == true
        PLOT.qualityplot = true;
    else
        PLOT.qualityplot = false;
    end
    if cfg.CollectFFT.correctpf == true
        PLOT.correctpf = true;
    else
        PLOT.correctpf = false;
    end
    if isfield(cfg.CollectFFT,'Limits')
        PLOT.Limits = cfg.CollectFFT.Limits;
    end
    PLOT.freqlabel = calc.freqlabel;
    
    % plot spectra
    if skipfig == 0
        set(f,'Visible','on','Units','pixels');
        pos = get(f,'position');
        shift = 70;
        U.Hprint = [];
        if Nplot < MaxPlot
            U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'pixels',...
                'position', [pos(3)-shift,2,70,25],...
                'string','next','fontsize',12,'Callback',@go_next,...
                'TooltipString','next');
        else
            U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'pixels',...
                'position', [pos(3)-shift,2,70,25],...
                'string','end','fontsize',12,'Callback',@go_end,...
                'TooltipString','end');
        end
        shift = shift + 75;
        if Nplot > 1
            U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'pixels',...
                'position', [pos(3)-shift,2,70,25],...
                'string','previous','fontsize',12,'Callback',@go_previous,...
                'TooltipString','previous');
            shift = shift + 75;
        end
        if skipedit == 0
            U.Hprint(1,end+1) = uicontrol('style', 'pushbutton', 'units', 'pixels',...
                'position', [pos(3)-shift,2,70,25],...
                'string','set NaN','fontsize',12,'Callback', ...
                @(hObj,evd)set_NaN(hObj,evd,cfg.CollectFFT),...
                'TooltipString','set NaN');
        end
        if isfield(PLOT.R,'good') & ~isempty(PLOT.R.good)
            GoodString = 'Good';
            if isfield(PLOT.R,'bchans')
                GoodString = [GoodString ' (bad chans: ' num2str(PLOT.R.bchans)]; %#ok<AGROW>
            end
            if isfield(PLOT.R,'badact')
                if ~strcmp(GoodString,'Good')
                    GoodString = [GoodString ' / bad activations: ' num2str(PLOT.R.badact(1,1))]; %#ok<AGROW>
                else
                    GoodString = [GoodString ' (bad activations: ' num2str(PLOT.R.badact(1,1))]; %#ok<AGROW>
                end
            end
            if ~strcmp(GoodString,'Good')
                GoodString = [GoodString ')']; %#ok<AGROW>
            end
            Bg = uicontrol('Parent',f,'Position',[5 2 270 25],'Style','checkbox', ...
                'HorizontalAlignment','left','BackgroundColor',[1 1 1],'fontsize',9, ...
                'String',GoodString,'Callback',@set_good);
            U.Hprint(1,end+1) = Bg;
            if PLOT.R.good == true
                set(Bg,'Value',1);
            else
                set(Bg,'Value',0);
            end
        else
            Bg = [];
        end
        if skipedit == 0
            set(f,'windowbuttondownfcn',@(hObj,evd)select_pf(hObj,evd,cfg.CollectFFT))  %Passing points to be edited
        end
        lab_plot_spectra(PLOT);
        set(f,'UserData',U);
        uiwait(f);
        if ~ishandle(f)
            flagexit = 1;
        end
    else
        lab_plot_spectra(PLOT);
    end
    
    if skipfig == 1
        lab_print_figure(fullfile(calc.resultplots_path,['BackgroundActivity_' PLOT.patient '.jpg']),f);
        Nplot = Nplot + 1;
        if Nplot > MaxPlot
            flagexit = 1;
        end
    elseif dowriting == true & flagexit == 0
        lab_print_figure(fullfile(calc.resultplots_path,['BackgroundActivity_' PLOT.patient '.jpg']),f);
    end
    clf;
    clearvars PLOT Bn Bp
end
if ishandle(f)
    close(f);
end
pause(0.2);

    function select_pf(H,H2,settings)
        ax =get(gcbf,'CurrentAxes');
        PF = get(ax,'currentpoint');
        if PF(1) > size(PLOT.T.SpectAll,2)
            return
        end
        if isfield(PLOT.T,'domedian')
            settings.domedian = PLOT.T.domedian;
        end
        if isfield(PLOT.T,'mappings')
            settings.mappings = PLOT.T.mappings;
        end
        if isfield(PLOT.T,'mappingBA')
            settings.mappingBA = PLOT.T.mappingBA;
        end
        if isfield(PLOT.T,'Valid')
            settings.Valid = PLOT.T.Valid;
        end
        if PF(1) >= 0
            settings.PF = PF(1);
            [R,T] = lab_collect_spectraldata_calcPFMF(PLOT.T.SpectAll,PLOT.T.SpectF,settings,PLOT.T);
            PLOT.R = R;
            PLOT.T = T;
            delete(ax);
            lab_plot_spectra(PLOT);
            if ~isempty(Bg) & PLOT.R.good == true
                set(Bg,'Value',1);
            else
                set(Bg,'Value',0);
            end
        end
    end
    
    function set_good(H,H2)
        if PLOT.R.good == true
            PLOT.R.good = false;
            set(H,'Value',0);
        else
            PLOT.R.good = true;
            set(H,'Value',1);
        end
    end

    function set_NaN(H,H2,settings)
        ax =get(gcbf,'CurrentAxes');
        if isfield(PLOT.T,'domedian')
            settings.domedian = PLOT.T.domedian;
        end
        if isfield(PLOT.T,'mappings')
            settings.mappings = PLOT.T.mappings;
        end
        if isfield(PLOT.T,'mappingBA')
            settings.mappingBA = PLOT.T.mappingBA;
        end
        if isfield(PLOT.T,'Valid')
            settings.Valid = PLOT.T.Valid;
        end
        settings.PF = NaN;
        [R,T] = lab_collect_spectraldata_calcPFMF(PLOT.T.SpectAll,PLOT.T.SpectF,settings,PLOT.T);
        PLOT.R = R;
        PLOT.T = T;
        delete(ax);
        lab_plot_spectra(PLOT);
        if ~isempty(Bg) & PLOT.R.good == true
            set(Bg,'Value',1);
        else
            set(Bg,'Value',0);
        end
    end
    
    function set_Limits(H,H2)
        if ~isfield(PLOT,'Limits') | ~isfield(PLOT.Limits,'y1')
            tmp = get(gca,'YLim');
            PLOT.Limits.y1 = tmp(1);
            PLOT.Limits.y2 = tmp(2);
            clearvars tmp
        end
        Prompt = {'Limit low','y1';'Limit high','y2'};
        Formats.type = 'edit';
        Formats.format = 'float';
        Formats.limits = [0 inf];
        PLOT.Limits = inputsdlg(Prompt,'Channels',[Formats;Formats],PLOT.Limits);
        set(gca,'YLim',[PLOT.Limits.y1 PLOT.Limits.y2]);
    end

    function go_next(H,H1)
        result.R(1,Nplot) = PLOT.R;
        result.T(1,Nplot) = PLOT.T;
        if Nplot < MaxPlot
            Nplot = Nplot+1;
        end
        uiresume(gcbf);
    end

    function go_select(H,H1)
        result.R(1,Nplot) = PLOT.R;
        result.T(1,Nplot) = PLOT.T;
        selection = listdlg('PromptString','Select','Name','Select Plot','SelectionMode','single', ...
            'ListString',result.patient,'CancelString','None','ListSize',[250 280]);
        if ~isempty(selection)
            Nplot = selection;
        end
        uiresume(gcbf);
    end

    function go_end(H,H1)
        result.R(1,Nplot) = PLOT.R;
        result.T(1,Nplot) = PLOT.T;
        flagexit = 1;
        uiresume(gcbf);
    end

    function go_previous(H,H1)
        result.R(1,Nplot) = PLOT.R;
        result.T(1,Nplot) = PLOT.T;
        if Nplot > 1
            Nplot = Nplot-1;
        end
        uiresume(gcbf);
    end
end

