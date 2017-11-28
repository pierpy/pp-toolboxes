% Create dialog or menu for 'TAPEEG-Main'
%
% lab_select_main(domenu)
%       domenu:  1            create menu in active figure
%                0 or empty   create dialog to select
%
% Written by F. Hatz 2012 Neurology Basel

function lab_select_main(domenu)

if ~exist('domenu','var')
    domenu = 0;
end

MenuName = 'Main';
SelList = {'Load data','lab_load_eeg','','off'; ...
    'Process data','lab_tapeeg','','off'; ...
    'Quit','quit','','off'};

lab_select_menu(MenuName,SelList,domenu);

return