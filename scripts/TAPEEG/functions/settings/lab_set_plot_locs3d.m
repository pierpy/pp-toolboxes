function [settings,skipprocessing] = lab_set_plot_locs3d(settings)

global Main_Path
    
if ~exist('settings','var') | isempty(settings) | ~isfield(settings,'locfile')
    settings.locfile = '';
    settings.dotscolor = [0 1 0];
    settings.dotscolorI = [1 0 0];
    settings.gapfactor = 4;
    settings.sizedots = 4;
    settings.plotlabels = true;
    settings.labelsize = 16;
    settings.labeldistance = 1;
    settings.mrifile = '';
    settings.facecolor = [0.7 0.7 0.7];
    settings.alpha = 0.5;
    settings.plotbrain = true;
    settings.Bfacecolor = [0.7 0.7 0.7];
    settings.Balpha = 1;
end

List_locs = {};
if ~isempty(Main_Path) & exist(fullfile(Main_Path,'electrodes'),'dir')
    List_locs = lab_search(fullfile(Main_Path,'electrodes'),{'*.els','*.xyz','*.sfp','*.elc','*.spi'},true,true,1);
end
List_locs = cat(1,{'Select File'},List_locs(:));
List_locs = cat(1,{settings.locfile},List_locs(:));

Prompt = {};
Formats = {};

Prompt(end+1,:) = {'LOC-file','locfile'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input';
Formats(end,1).items = List_locs;
Formats(end,1).size = 250;
Formats(end,1).callback = {@load_locs,{'LOCS','LOC_file'},'locfile'};
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Color','dotscolor'};
Formats(end+1,1).type = 'color';
Formats(end,1).format = 'color';

Prompt(end+1,:) = {'Size', 'sizedots'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Plot labels','plotlabels'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Size', 'labelsize'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;

Prompt(end+1,:) = {'Distance', 'labeldistance'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf];
Formats(end,1).size = 40;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {' ', ''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'MRI-file','mrifile'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.hdr','MRI-file (*.hdr)'};
Formats(end,1).limits = [0 1];
Formats(end,1).size = 300;
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Color','facecolor'};
Formats(end+1,1).type = 'color';
Formats(end,1).format = 'color';

Prompt(end+1,:) = {'Alpha','alpha'};
Formats(end+1,1).type = 'range';
Formats(end,1).limits = [0 1];

Prompt(end+1,:) = {' ', ''};
Formats(end+1,1).type = 'text';
Formats(end,1).span = [1 5];

Prompt(end+1,:) = {'Plot Brain','plotbrain'};
Formats(end+1,1).type = 'check';
Formats(end,1).span = [1 3];

Prompt(end+1,:) = {'Color','Bfacecolor'};
Formats(end+1,1).type = 'color';
Formats(end,1).format = 'color';

Prompt(end+1,:) = {'Alpha','Balpha'};
Formats(end+1,1).type = 'range';
Formats(end,1).limits = [0 1];

[settings,skipprocessing] = inputsdlg(Prompt,'Plot electrodes 3d',Formats,settings);
pause(0.2);

    function LOC_file = load_locs(LOC_file)
        if ~exist(LOC_file,'file') | strcmp(LOC_file,'Select File')
            [LOC_file,LOC_filepath] = uigetfile('*.els;*.sfp;*.xyz','Select LOC-file');
            if LOC_file == 0
                return
            end
            LOC_file = fullfile(LOC_filepath,LOC_file);
        end
    end

end