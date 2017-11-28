% Function to save figures as image/pdf
% (wrapper for export_fig from mathworks file exchange)
%
% lab_print_figure(Filename,Fig,FFormat,Mode)
%
% Filename  = Filepath&Filename to store image
% Fig       = Figure handle (optional)
% Format    = File format (optional, if empty defined by extension of Filename)
% Mode      = simple (simple print command without using export_fig)
%             advanced (export_fig)
%
% written by F. Hatz 2014, University Hospital Basel

function lab_print_figure(Filename,Fig,Format,Mode,settings)

if ~exist('settings','var')
    settings = [];
end
if exist('Filename','var') & ~isempty(Filename)
    if ~exist('Mode','var')
        Mode = 'simple';
    end
    if exist('Format','var') & ~isempty(Format)
        [~,Filepath,~,FilenameS] = lab_filename(Filename);
        Filename = [FilenameS '.' lower(Format)];
    else
        [Filename,Filepath,Format,FilenameS] = lab_filename(Filename);
    end
    Filename = regexprep(Filename,{'::',':',';','/','\',','},'_');
    if isempty(Filepath)
        Filepath = pwd;
    end
else
    if ~exist('Mode','var')
        Mode = 'advanced';
    end
    settings = lab_set_print_figure(settings);
    if isempty(settings)
        return
    end
    Filename = settings.Filename;
    [Filename,Filepath,Format,FilenameS] = lab_filename(Filename);
end
if ~exist('Fig','var') | isempty(Fig)
    Fig = gcf;
end

% Goto folder
cd(Filepath);

% Find Layout-Objects to set invisible
U = get(Fig,'UserData');
if isfield(U,'Hprint') & ~isempty(U.Hprint)
    for i = 1:length(U.Hprint)
        if ishandle(U.Hprint(i))
            set(U.Hprint(i),'Visible','off');
        end
    end
end

if exist(fullfile(Filepath,Filename),'file') & (~isfield(settings,'addpdf') | settings.addpdf == false)
    delete(fullfile(Filepath,Filename));
end

if strcmp(Mode,'simple')
    switch Format
        case 'tif'
            dformat = '-dtiff';
        case 'tiff'
            dformat = '-dtiff';
        case 'bmp'
            dformat = '-dbmp';
        case 'emf'
            dformat = '-dmeta';
        case 'pdf'
            dformat = '-dpdf';
        case 'eps'
            dformat = '-depsc';
        case 'png'
            dformat = '-dpng';
        otherwise
            dformat = '-djpeg';
            Filename = [FilenameS '.jpg'];
    end
    warning off %#ok<WNOFF>
    print(Fig,dformat,fullfile(Filepath,Filename));
    warning on %#ok<WNON>
else
    set(Fig, 'Units', 'centimeters','PaperUnits','centimeters')
    pos = get(Fig,'Position');
    set(Fig,'PaperSize', [pos(3) pos(4)]);
    set(Fig,'PaperPositionMode', 'manual');
    set(Fig,'PaperPosition', [0 0 pos(3) pos(4)]);
    
    if ~isfield(settings,'antialiasing')
        settings.antialiasing = 4;
    end
    if ~isfield(settings,'quality')
        settings.quality = 80;
    end
    if ~isfield(settings,'resolution')
        settings.resolution = 300;
    end
    if ~isfield(settings,'renderer')
        settings.renderer = lower(get(gcf,'Renderer'));
        if ~strcmp(settings.renderer,'opengl')
            settings.renderer = 'painters';
        end
    end
    
    switch Format
        case 'tif'
            dformat = '-tif';
        case 'tiff'
            dformat = '-tif';
        case 'bmp'
            dformat = '-bmp';
        case 'emf'
            dformat = '-emf';
        case 'pdf'
            dformat = '-pdf';
        case 'eps'
            dformat = '-eps';
        case 'png'
            dformat = '-png';
        otherwise
            dformat = '-jpg';
            Filename = [FilenameS '.jpg'];
    end
    warning off %#ok<WNOFF>
    if strcmp(dformat,'-emf')
        print(Fig,'-dmeta',fullfile(Filepath,Filename));
    else
        parameters = {['-a' num2str(settings.antialiasing)], ...
            ['-r' num2str(settings.resolution)],['-q' num2str(settings.quality)]};
        if ~isfield(settings,'docrop') | settings.docrop == false
            parameters{1,end+1} = '-nocrop';
        end
        if ~strcmp(settings.renderer,'auto')
            parameters{1,end+1} = ['-' settings.renderer];
        end
        if isfield(settings,'transparent')
            parameters{1,end+1} = '-transparent';
        end
        if isfield(settings,'native')
            parameters{1,end+1} = '-native';
        end
        if isfield(settings,'colorspace')
            parameters{1,end+1} = ['-' lower(settings.colorspace)];
        end
        if isfield(settings,'addpdf') & settings.addpdf == true
            parameters{1,end+1} = '-append';
        end
        export_fig(Fig,fullfile(Filepath,Filename),dformat,parameters{:});
    end
    warning on %#ok<WNON>
end

if isfield(U,'Hprint') & ~isempty(U.Hprint)
    for i = 1:length(U.Hprint)
        if ishandle(U.Hprint(i))
            set(U.Hprint(i),'Visible','on');
        end
    end
end

return