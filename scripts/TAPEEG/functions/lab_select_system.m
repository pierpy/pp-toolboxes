% Create dialog or menu for 'System'
%
% lab_select_settings(domenu)
%       domenu:  1            create menu in active figure
%                0 or empty   create dialog to select
%
% Written by F. Hatz 2012 Neurology Basel

function lab_select_system(domenu)

if ~exist('domenu','var')
    domenu = 0;
end

MenuName = 'System';
SelList = {'Set System-Font',{'lab_selectfont_fixedwidth',1},'','off'; ...
    'Set Debug-Mode',{'lab_debug'},'','off'; ...
    'Set lag for amplifier','lab_define_egi_lag','','off'};

lab_select_menu(MenuName,SelList,domenu);

return