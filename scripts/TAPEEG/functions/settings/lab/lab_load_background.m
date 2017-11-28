% written by F. Hatz 2013

function BA = lab_load_background(BA,list,locs)

if ~exist('locs','var')
    locs = [];
end
if ~exist('list','var')
    list = [];
end
if ~exist('BA','var')
    BA = [];
end
if isempty(list) & isfield(locs,'labels')
    list = locs.labels';
end

if isempty(BA) & ~isempty(list)
    BA = lab_get_bactivity(size(list,1));
elseif isempty(BA)
    BA = [];
else
    BA = BA(:)';
end

settings.BA = BA;
if exist('BA_file','var')
    settings.BA_file = BA_file;
else
    settings.BA_file = [];
end
Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Background Channels - File','BA_file'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.xls;*.xlsx','Background-File (.xls,.xlsx)'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).callback =  {@read_BA,'BA','BA_file','BA',size(list,1)};
Formats(end,1).size = 270;

Prompt(end+1,:) = {'Channels','BA'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).size = 360;
Formats(end,1).limits = [-inf inf];
if ~isempty(locs)
    Formats(end,1).enable = 'inactive';
    Formats(end,1).callback =  {@select_BA,'BA','BA',locs};
elseif ~isempty(list)
    Formats(end,1).items = list;
end

if ~isempty(list)
    Prompt(end+1,:) = {'Save',''};
    Formats(end+1,1).type = 'button';
    Formats(end,1).style = 'pushbutton';
    Formats(end,1).size = [100 25];
    Formats(end,1).callback = {@save_BA,'BA_file','BA',size(list,1)};
end

[settings,Cancelled] = inputsdlg(Prompt,'Background Channels',Formats,settings);
if Cancelled == 1
    if ~isempty(list)
        BA = 1:size(list,1);
    else
        BA = [];
    end
    return
end
BA = settings.BA;

   function BA = read_BA(BA_file,BA,numchans)
        if ~isempty(BA_file) & exist(BA_file,'file')
            if ispc
                [BA,~,xlsinput] = xlsread(BA_file,1);
            else
                [BA,~,xlsinput] = xlsread(BA_file,1,'','basic');
            end
            if size(xlsinput,2) == 2 & numchans == 0
                List = cell(size(xlsinput,1),1);
                for i = 1:size(List,1)
                    List{i,1} = num2str(xlsinput{i,1});
                end
                BA = listdlg('PromptString','Number of channels','SelectionMode','single','ListString',List');
                if ~isempty(BA)
                    BA = str2num(xlsinput{BA,2}); %#ok<ST2NM>
                    if ischar(BA)
                        BA = str2num(BA); %#ok<ST2NM>
                    end
                else
                    BA = [];
                end
            elseif size(xlsinput,2) == 2 & ~isempty(find(BA==numchans,1))
                BA = xlsinput{find(BA==numchans,1),2};
                if ischar(BA)
                    BA = str2num(BA); %#ok<ST2NM>
                end
            elseif size(xlsinput,2) == 2 & ischar(xlsinput{1,2})
                BA= [];
            elseif size(BA,1) > 1 & size(BA,2) == 1
                BA = BA';
            elseif size(BA,1) > 1
                BA = BA(1,:);
            end
       end
   end
   
   function BA_file = save_BA(BA,numchans)
       if isempty(BA)
           return
       end
       [BA_file,BA_filepath] = uiputfile('*.xls;*.xlsx','Select File to store');
       if BA_file ~= 0
           BA_file = fullfile(BA_filepath,BA_file);
           if numchans > 0
               xlsout{1,1} = numchans;
               xlsout{1,2} = num2str(BA);
           else
               xlsout = BA;
           end
           lab_write_xls(BA_file,xlsout);
       else
           BA_file = '';
       end
   end

   function BA = select_BA(BA,locs)
       if isfield(locs,'x')
           BA = BA(BA<=length(locs.x));
       end
       plot.indexed = BA;
       plot.Color = [1 1 1];
       plot.ColorIdx = [1 0 0];
       plot.LOCS = locs;
       plot.Title = 'Background channels';
       BA = lab_plot_locs(plot,1,0,0);
   end
end