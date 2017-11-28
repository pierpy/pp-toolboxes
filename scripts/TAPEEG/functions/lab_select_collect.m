% Create dialog or menu for 'Collect'
%
% lab_select_collect(domenu)
%       domenu:  1            create menu in active figure
%                0 or empty   create dialog to select
%
% Written by F. Hatz 2012 Neurology Basel

function lab_select_collect(domenu)

if ~exist('domenu','var')
    domenu = 0;
end

MenuName = 'Collect';
SelList = {'Bad channels','lab_collect_badchannels',''; ...
    'Frequency data','lab_collect_spectraldata',''; ...
    'Connectivity data','lab_collect_connectivity',''; ...
    'Graph analysis','lab_collect_graphanalysis',''; ...
    'Microstates','lab_collect_microstates',''; ...
    'Distance matrices','lab_collect_distances',''; ...
    'AVG Sweeps','lab_collect_avgsweeps',''};

lab_select_menu(MenuName,SelList,domenu);

return