function lab_topos2video

[Filename,Filepath] = uigetfile('.sef');
if isempty(Filename) | isnumeric(Filename)
    return
end
[~,~,~,FilenameS] = lab_filename(Filename);

% create output folder
Video_file = FilenameS;
warning off %#ok<WNOFF>
Video_filepath = fullfile(Filepath,['TopoVideo_' FilenameS]);
mkdir (Video_filepath);
warning on %#ok<WNON>

[data,header] = lab_read_data(fullfile(Filepath,Filename));
if ~isfield(header,'locs') | isempty(header.locs)
    disp('Abort, no LOCS information')
    return
end

settings.LOCS = header.locs;
if isfield(header,'numdatachannels')
    [data,header] = lab_reduce_channels(data,header,1:header.numdatachannels); %#ok<NASGU>
end

if size(data,1) ~= size(settings.LOCS.x,2)
    disp('Abort, wrong LOCS information')
    return
end

disp('   Write video')
settings.PLOT_file = [];
settings.close = 0;
settings.AddPlot = 0;
settings.startframe = 1;
settings.stopframe = size(data,2);
settings.NoMenu = true;
range = max(abs(data(:)));
PLOT.MinValue = -range;
PLOT.MaxValue = range;
PLOT.Color = lab_create_cmap('bluered');
Prompt = {'Start Frame', 'startframe';'Stop Frame', 'stopframe'};
Formats.type = 'edit';
Formats.format = 'integer';
Formats.limits = [1 settings.stopframe-1];
Formats.size = 80;
Formats(2).type = 'edit';
Formats(2).format = 'integer';
Formats(2).limits = [2 settings.stopframe];
Formats(2).size = 80;
[settings,Cancelled] = inputsdlg(Prompt,'Video Length',Formats,settings);
pause(0.2);
if Cancelled == 1
    return
end
data = data(:,settings.startframe:settings.stopframe);
settings.nframe = size(data,2);
Mov = VideoWriter(fullfile(Video_filepath,Video_file),'MPEG-4');
open(Mov);
tmp.colormap = [];
progressbar('store video')
settings.handleF = figure('Visible','off','Color',[1 1 1],'Renderer','zbuffer','Menubar','none');
for i = 1:settings.nframe
    settings = lab_plot_chans(data(:,i),PLOT,settings);
    tmp.cdata = zbuffer_cdata(settings.handleF);
    writeVideo(Mov,tmp);
    delete(settings.handleL);
    progressbar(i/settings.nframe);
end

end

function cdata = zbuffer_cdata(hfig)
   % Get CDATA from hardcopy using zbuffer
   % Need to have PaperPositionMode be auto
   orig_mode = get(hfig, 'PaperPositionMode');
   set(hfig, 'PaperPositionMode', 'auto');
   cdata = hardcopy(hfig, '-Dzbuffer', '-r0');
   % Restore figure to original state
   set(hfig, 'PaperPositionMode', orig_mode);
end

