% Create dialog or menu for 'Settings'
%
% lab_select_settings(domenu)
%       domenu:  1            create menu in active figure
%                0 or empty   create dialog to select
%
% Written by F. Hatz 2012 Neurology Basel

function lab_select_settings(domenu)

global Main_Path

if ~exist('domenu','var')
    domenu = 0;
end

MenuName = 'Settings';
SelList = {'Show settings (save disabled)',{'lab_edit_settings',[],'Read'},'','off'; ...
    'Edit settings',{'lab_edit_settings',[],'Edit'},'','off'; ...
    'EEG/MEG processing',{'lab_edit_settings',[],'Create'},'Create new settings','off'; ...
    'MRI processing',{'lab_edit_settings',[],'CreateMRI'},'Create new settings','off'; ...
    'Matrix processing',{'lab_edit_settings',[],'CreateMatrix'},'Create new settings','off'; ...
    'Collect Spectra',{'lab_edit_settings',[],'CreateCollectSpectra'},'Create new settings','off'; ...
    'Collect Graph',{'lab_edit_settings',[],'CreateCollectGraph'},'Create new settings','off'; ...
    'Collect Connectivity',{'lab_edit_settings',[],'CreateCollectConnectivity'},'Create new settings','off'; ...
    'Collect Microstates',{'lab_edit_settings',[],'CreateCollectMicrostates'},'Create new settings','off'; ...
    'Edit/Create auto-settings',{'lab_edit_settings',[],Main_Path},'','on'};

lab_select_menu(MenuName,SelList,domenu);

return