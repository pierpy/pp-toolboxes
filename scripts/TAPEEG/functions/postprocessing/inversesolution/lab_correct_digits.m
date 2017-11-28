% Delete outliers in digits-information. Input usually is fiff-file
%
% function digits = lab_correct_digits(digits)
%
% written by F. Hatz 2012

function digits = lab_correct_digits(digits)

if ~exist('digits','var')
    [LOC_file,LOC_filepath]=uigetfile('*.xyz;*.els;*.spi;*.fif;*.fiff;*.mff','Select digits-file');
    if LOC_file == 0
        digits = [];
        return
    else
        Filelist = {fullfile(LOC_filepath,LOC_file)};
        clearvars LOC_file LOC_filepath
    end
elseif iscell(digits) & ~isempty(digits) & ischar(digits{1,1})
    Filelist = digits;
elseif isnumeric(digits) & size(digits,2) == 3
    Filelist{1} = digits;
elseif isstruct(digits) & isfield(digits,'digits')
    Filelist{1} = digits;
else
    digits = [];
    return
end

for filenr = 1:size(Filelist,2)
    if isnumeric(Filelist{1,filenr})
        digits = Filelist{1,filenr};
    elseif isstruct(Filelist{1,filenr}) & isfield(Filelist{1,filenr},'digits')
        LOCS = Filelist{1,filenr};
        digits = LOCS.digits;
    elseif ischar(Filelist{1,filenr}) & exist(Filelist{1,filenr},'file')
        [~,LOC_filepath,LOC_format,LOC_fileS] = lab_filename(Filelist{1,filenr});
        if ~exist('LOC_format','var') | ~strcmp(LOC_format','xyz')
            LOC_format = 'els';
        end
        LOCS = lab_read_locs(Filelist{1,filenr});
        if ~isempty(LOCS) & isfield(LOCS,'digits') & ~isempty(LOCS.digits)
            digits = LOCS.digits;
        else
            digits = [];
        end
    else
        digits = [];
    end
    
    if ~isempty(digits)
        digitsO = digits;
        
        x = digits(:,1);
        y = digits(:,2);
        z = digits(:,3);
        
        f = figure('Color',[1 1 1],'InvertHardCopy','off','NumberTitle','off','Menubar','none');
        a = axes;
        set(a,'Visible','off','Parent',f,...
            'Position',[0.01 0.01 0.99 0.99],'xtick',[],'ytick',[], ...
            'PlotBoxAspectRatio',[1 1 1],...
            'DataAspectRatio',[1 1 1]);
        
        %Plotting for the points to be edited
        hold on
        pH = plot3(a,x,y,z,'ro');
        %Turning it off in case of a scatterplot
        hold off
        set(f,'windowbuttondownfcn',@select_point)  %Passing points to be edited
        
        pos = get(f,'position');
        btnOK = uicontrol('style', 'pushbutton', 'units', 'pixels',...
            'position', [pos(3)-50,2,50,30],...
            'string','ok','fontsize',14,'Callback','uiresume(gcbf);',...
            'TooltipString','ok');
        btnCancel = uicontrol('style', 'pushbutton', 'units', 'pixels',...
            'position', [pos(3)-145,2,90,30],...
            'string','Cancel','fontsize',14,'Callback','close(gcbf);',...
            'TooltipString','Cancel');
        btnRot = uicontrol('style', 'checkbox', 'units', 'pixels',...
            'position', [2,2,80,25],...
            'string','Rotate','fontsize',10,'Callback','rotate3d',...
            'TooltipString','rotate on/off');
        uiwait(f)
        
        if ishandle(f)
            close(f);
        else
            digits = digitsO;
        end
        if exist('LOC_filepath','var')
            LOCS.digits = digits;
            lab_write_locs(fullfile(LOC_filepath,[LOC_fileS '.' LOC_format]),LOCS);
            clearvars LOC_fileS LOC_filepath LOC_format LOCS
        elseif exist('LOCS','var')
            LOCS.digits = digits;
            digits = LOCS;
            clearvars LOCS
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
    function select_point(varargin)
        [~,~,vi] = select3d;
        xp = get(pH,'Xdata'); %Getiing coordinates of points to be edited
        yp = get(pH,'Ydata');
        zp = get(pH,'Zdata');
        xp = xp(1,setdiff(1:size(xp,2),vi));
        yp = yp(1,setdiff(1:size(yp,2),vi));
        zp = zp(1,setdiff(1:size(zp,2),vi));
        set(pH,'Xdata',xp);
        set(pH,'Ydata',yp);
        set(pH,'Zdata',zp);       
        digits = digits(setdiff(1:size(digits,1),vi),:);
        drawnow
    end
end
