% Dialog to select orientation of input mri, used by lab_plot_orthoslides
%
% written by F. Hatz

function [orient,Cancelled] = lab_get_orient(orient)

if ~exist('orient','var') | isempty(orient)
    orient = [1 2 3];
end
settings.x = orient(1);
settings.y = orient(2);
settings.z = orient(3);
choice = {'From Left to Right','From Posterior to Anterior', ...
    'From Inferior to Superior','From Right to Left', ...
    'From Anterior to Posterior','From Superior to Inferior',''};

Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Orientation of the original image:',''};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'X Axes','x'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = choice;
Formats(end,1).callback = {@check_orient,'@ALL','@ALL',1};

Prompt(end+1,:) = {'Y Axes','y'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = choice;
Formats(end,1).callback = {@check_orient,'@ALL','@ALL',2};

Prompt(end+1,:) = {'Z Axes','z'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = choice;
Formats(end,1).callback = {@check_orient,'@ALL','@ALL',3};

Options.WindowStyle = 'normal';
Options.ForceWindowStyle = true;
[settings,Cancelled] = inputsdlg(Prompt,'Set original orientation',Formats,settings,Options);
if Cancelled == 0
   orient = [settings.x settings.y settings.z];
   if ~isempty(find(orient==7,1))
       orient = [1 2 3];
   end
end

end

function settings = check_orient(settings,number)

orig = {settings.x,settings.y,settings.z};
tmp = setdiff(1:3,number);
new = [orig{number} 0];
if new(1) > 3 & new(1) < 7
    new(2) = new(1) - 3;
elseif new(1) < 7
    new(2) = new(1) + 3;
else
    new(1) = 0;
end
if orig{tmp(1)} == new(1) | orig{tmp(1)} == new(2)
    orig{tmp(1)} = 7;
end
if orig{tmp(2)} == new(1) | orig{tmp(2)} == new(2)
    orig{tmp(2)} = 7;
end
settings.x = orig{1};
settings.y = orig{2};
settings.z = orig{3};

end