% Create dialog or menu for TAPEEG
%
% Written by F. Hatz 2012 Neurology Basel

function lab_select_menu(MenuName,SelList,domenu)

if exist('domenu','var') & domenu == 1
    m1 = uimenu(gcf,'Label',MenuName);
    tmp = {};
    tmp2 = {};
    for i = 1:size(SelList,1)
        if size(SelList,2) > 2 & ~isempty(SelList{i,3})
            index = find(strcmp(tmp,SelList{i,3}));
            if isempty(index)
                tmp2{end+1} = uimenu(m1,'Label',SelList{i,3}); %#ok<AGROW>
                m2 = tmp2{end};
                tmp{end+1} = SelList{i,3}; %#ok<AGROW>
            else
                m2 = tmp2{index};
            end
            clearvars index
            if size(SelList,2) > 3 & strcmp(SelList{i,4},'on')
                uimenu(m2,'Label',SelList{i,1},'Callback', ...
                    set_command(SelList{i,2}),'Separator','on');
            else
                uimenu(m2,'Label',SelList{i,1},'Callback', ...
                    set_command(SelList{i,2}));
            end
        else
            if size(SelList,2) > 3 & strcmp(SelList{i,4},'on')
                uimenu(m1,'Label',SelList{i,1},'Callback', ...
                    set_command(SelList{i,2}),'Separator','on');
            else
                uimenu(m1,'Label',SelList{i,1},'Callback', ...
                    set_command(SelList{i,2}));
            end
        end
    end
    clearvars tmp tmp2
else
    for i = 1:size(SelList,1)
        if size(SelList,2) > 2 & ~isempty(SelList{i,3})
            strlist{i,1} = [SelList{i,3} ' - ' SelList{i,1}]; %#ok<AGROW>
        else
            strlist{i,1} = SelList{i,1}; %#ok<AGROW>
        end
    end
    [selection] = listdlg('PromptString','Select','Name',MenuName,'SelectionMode','single', ...
        'ListString',[{'none'} strlist'],'CancelString','None','ListSize',[250 280]);
    if isempty(selection) | selection == 1
        return
    else
        eval([set_command2(SelList{selection-1,2}) ';']);
        clearvars selection SelList
        disp('Finished')
    end
end

end

function command = set_command(command)
   if ischar(command) & ~isempty(command)
       if exist(command) %#ok<EXIST>
           command = ['close;pause(0.2);' command ';lab_show_start;'];
       else
           command = 'msgbox(''not released yet'');';
       end
   elseif iscell(command) & ~isempty(command)
       if exist(command{1}) %#ok<EXIST>
           if length(command) > 1
               tmp = '(';
               for i = 2:length(command)
                   if isempty(command{i})
                       tmp = [tmp '[]']; %#ok<AGROW>
                   elseif isnumeric(command{i})
                       tmp = [tmp num2str(command{i})]; %#ok<AGROW>
                   elseif ischar(command{i})
                       tmp = [tmp '''' command{i} '''']; %#ok<AGROW>
                   end
                   if i == length(command)
                       tmp = [tmp ')']; %#ok<AGROW>
                   else
                       tmp = [tmp ',']; %#ok<AGROW>
                   end
               end
               command = [command{1} tmp];
           else
               command = command{1};
           end
           command = ['close;pause(0.2);' command ';lab_show_start;'];
       else
           command = 'msgbox(''not released yet'');';
       end
   else
       command = 'msgbox(''internal error'');';
   end
end

function command = set_command2(command)
   if ischar(command) & ~isempty(command)
       if exist(command) %#ok<EXIST>
           command = [command ';'];
       else
           command = 'msgbox(''not released yet'');';
       end
   elseif iscell(command) & ~isempty(command)
       if exist(command{1}) %#ok<EXIST>
           if length(command) > 1
               tmp = '(';
               for i = 2:length(command)
                   if isempty(command{i})
                       tmp = [tmp '[]']; %#ok<AGROW>
                   elseif isnumeric(command{i})
                       tmp = [tmp num2str(command{i})]; %#ok<AGROW>
                   elseif ischar(command{i})
                       tmp = [tmp '''' command{i} '''']; %#ok<AGROW>
                   end
                   if i == length(command)
                       tmp = [tmp ')']; %#ok<AGROW>
                   else
                       tmp = [tmp ',']; %#ok<AGROW>
                   end
               end
               command = [command{1} tmp];
           else
               command = command{1};
           end
           command = [command ';'];
       else
           command = 'msgbox(''not implemented yet'');';
       end
   else
       command = 'msgbox(''internal error'');';
   end
end