% settings for detect bad
%
% lab_debug
%
% F. Hatz 2013

function lab_debug

global dodebug

if ~isfield(dodebug,'bad_container')
    dodebug.bad_container = false;
    dodebug.bad_container_folder = '';
end

Prompt = cell(0,2);
Formats = {};

Prompt(end+1,:) = {'Write Container-Files for Bad-Detection','bad_container'};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Store-Folder','bad_container_folder'};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = 300;

[dodebug,Cancelled] = inputsdlg(Prompt,'Debug Mode',Formats,dodebug);
if isempty(dodebug) | Cancelled == 1
    dodebug = [];
    return
end
pause(0.2);

end