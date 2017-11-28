%INPUTSDLG DEMO (Enhanced input dialog box with multiple data types)

% Written by: Takeshi Ikuma
% Last Updated: May 5 2010
%
% Updated for additional functions F. Hatz 2013

Title = 'INPUTSDLG Demo Dialog';

%%%% SETTING DIALOG OPTIONS
Options = [];
Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'on';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 4; % Horizontal dimension in fields

Prompt = {};
Formats = {};

Prompt(end+1,:) = {['This demo illustrates every type of control that can be placed by INPUTSDLG function ' ... 
   'and demonstrates how Formats input can be used to layout these controls.'],[],[]};
Formats(end+1,1).type = 'text';
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 4]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'Bidder''s Name', 'Name',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).size = 200; % automatically assign the height

Prompt(end+1,:) = {'Bidder''s SSN (no space or hyphen)', 'SSN',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'integer';
Formats(end,1).limits = [0 999999999]; % 9-digits (positive #)
Formats(end,1).size = 80;

Prompt(end+1,:) = {'Bidding Price', 'Price',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).limits = [0 inf]; % non-negative decimal number

Prompt(end+1,:) = {'Enable enhanced mode' 'EnableEnhancedMode',[]};
Formats(end+1,1).type = 'check';

Prompt(end+1,:) = {'Bidder''s Bio File','BioFile',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).items = {'*.bio','Biography File (*.bio)';'*.*','All Files'};
Formats(end,1).limits = [0 1]; % single file get
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 3];  % item is 1 field x 3 fields

Prompt(end+1,:) = {'Action','Action',[]};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'togglebutton';
Formats(end,1).items = {'Bid';'Decline';'Pass'};
Formats(end,1).span = [3 1];  % item is 3 fields x 1 field

Prompt(end+1,:) = {'Bidder''s Data Folder','DataFolder',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'dir';
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 3];  % item is 1 field x 3 fields

Prompt(end+1,:) = {'Save Bidding History To','SaveFile',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).limits = [1 0]; % use uiputfile
Formats(end,1).size = [-1 0];
Formats(end,1).span = [1 3];  % item is 1 field x 3 fields

Prompt(end+1,:) = {'Select Item Files','ItemFiles',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'file';
Formats(end,1).limits = [0 5]; % multi-select files
Formats(end,1).size = [-1 -1];
Formats(end,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
Formats(end,1).span = [1 3];  % item is 1 field x 3 fields

% Prompt(end+1,:) = {'Image:','Image',[]};
% Formats(end+1,1).type = 'image';
% Formats(end,1).size = [200 150];
% Formats(end,1).items = {'path to image file'}; % put here path to image-file (only used if no image in default answer)
% Formats(end,1).span = [3 1];

Prompt(end+1,:) = {'Image:','Image',[]};
Formats(end+1,1).type = 'image';
Formats(end,1).format = 'matrix';
Formats(end,1).size = [200 150];
Formats(end,1).items = {repmat((0:1/255:1)',1,3)}; % put here path to image-file (only used if no image in default answer)
Formats(end,1).span = [3 1];

Prompt(end+1,:) = {'Choose Currency','MoneyUnit',[]};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'radiobutton';
Formats(end,1).items = {'U.S. Dollar' 'Euro';'Japanese Yen' ''};
Formats(end,1).span = [2 1];  % item is 2 field x 1 fields

Prompt(end+1,:) = {'Item Table','Table',[]};
Formats(end+1,1).type = 'table';
Formats(end,1).format = 'table';
Formats(end,1).items = {{'value','left','right','select'},{'Row1','Row2'}};
Formats(end,1).size = [373 73];
Formats(end,1).span = [2 2];  % item is 2 field x 1 fields

Prompt(end+1,:) = {'Bidding Rate','BidRate',[]};
Formats(end+1,1).type = 'range';
Formats(end,1).limits = [0 1];

Prompt(end+1,:) = {'Result','Result',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'result';
Formats(end,1).callback = {@inputsdlg_demo_table, ... % string with script-name or function_handle 
                        [], ... & output variable (empty = own variable)
                        'Result', ... % take own variable (Answer.Result) as first parameter
                        {'value','left','right','select'}, ... & 2. parameter
                        'EditTable'}; % 3. parameter (=string)

Prompt(end+1,:) = {'Item List','List',[]};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'Black','White','Red','Blue','Green','Yellow','Orange'};

Prompt(end+1,:) = {'Item Color','Color',[]};
Formats(end+1,1).type = 'color';

Prompt(end+1,:) = {'Item Vector','Vector',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf]; % if not [0 1] data = numeric


Prompt(end+1,:) = {'Date (dd-mmm-yyyy)', 'Date',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'date';
Formats(end,1).limits = 'dd-mmm-yyyy'; % with time: 'dd-mmm-yyyy HH:MM:SS'
Formats(end,1).size = 100;

Prompt(end+1,:) = {'Item List (rank)','ListRank',[]};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).items = {'Black','White','Red','Blue','Green','Yellow','Orange'};

Formats(end+1,1).type = 'none'; % leave next place in dialog-array empty

Prompt(end+1,:) = {'Item Vector (strings)','VectorStrings',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [0 1]; % [0 1] data = strings

Prompt(end+1,:) = {'Memo','Memo',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'text';
Formats(end,1).limits = [0 9]; % default: show 9 lines
Formats(end,1).size = [-1 0];
Formats(end,1).span = [4 3];  % item is 4 fields x 3 fields

Prompt(end+1,:) = {'Item ListVector','ListVector',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).limits = [-inf inf]; % if not [0 1] data must be numeric otherwise string
Formats(end,1).items = {'val1','val2','val3','val4','val5','val6'}; % For edit list is shown with items

Prompt(end+1,:) = {'Item ListVector (strings)','ListVectorStrings',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'vector';
Formats(end,1).items = {'val1','val2','val3','val4','val5','val6'}; % For edit list is shown with items

Prompt(end+1,:) = {'','',[]};
Formats(end+1,1).type = 'text';

Prompt(end+1,:) = {'Auction Sites:','Site',[]};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'www.auction1.com';'www.auction2.com';'www.bidme.com';'www.bestvalu.com'};
Formats(end,1).limits = [0 4]; % multi-select
Formats(end,1).size = [140 80];
Formats(end,1).span = [3 1];  % item is 2 fields x 1 fields

Prompt(end+1,:) = {'Auction Sites (rank):','SiteRank',[]};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'listbox';
Formats(end,1).items = {'www.auction1.com','www.auction2.com','www.bidme.com','www.bestvalu.com'};
Formats(end,1).limits = [0 2]; % multi-select
Formats(end,1).size = [140 80];
Formats(end,1).span = [3 1];  % item is 2 fields x 1 fields

Prompt(end+1,:) = {'             X', 'X',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).size = 80;

Prompt(end+1,:) = {'X*Y','Button',[]};
Formats(end+1,1).type = 'button';
Formats(end,1).style = 'pushbutton';
Formats(end,1).size = [100 25];
Formats(end,1).callback = {@inputsdlg_demo_function, ... % function_handle
                        'Z', ...
                        'X', ...
                        'Y'};

Prompt(end+1,:) = {'             Y', 'Y',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).size = 80;

Prompt(end+1,:) = {'= ', 'Z',[]};
Formats(end+1,1).type = 'edit';
Formats(end,1).format = 'float';
Formats(end,1).size = 120;

Prompt(end+1,:) = {'', 'hidden',[]};
Formats(end+1,1).type = 'hidden';

%%%% SETTING DEFAULT STRUCT
DefAns.Name = 'John Smith';
DefAns.SSN = 123456789;
DefAns.Price = 99.99;
DefAns.EnableEnhancedMode = true;
d = dir;
files = strcat([pwd filesep],{d(~[d.isdir]).name});
DefAns.BioFile = files{1};
DefAns.Action = 2; % = 'Decline'
DefAns.DataFolder = pwd;
DefAns.SaveFile = files{2};
DefAns.ItemFiles = files(3:end);
DefAns.MoneyUnit = 3; % yen
DefAns.BidRate = 0.75;
DefAns.Table = {'val1',1,10,false;'val2',0,16,true};
DefAns.Memo = 'Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque elementum, dui sed sagittis vulputate, nulla tellus dapibus velit, pretium molestie odio lorem eget nisl. Cras volutpat gravida neque, vitae dictum eros aliquet non. Nam in sem sapien, sed condimentum justo. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Fusce ut dignissim elit. Aenean ac massa arcu. Pellentesque rhoncus fermentum tortor et auctor. Proin eu velit dolor, eu semper nunc. Praesent in lacinia orci. Etiam ac diam nibh. Ut aliquet sapien sed metus blandit vitae aliquam enim ullamcorper. In nec ligula id quam interdum porttitor. Sed tempus turpis ut est convallis et interdum lectus dapibus. In dictum pulvinar erat id ornare.';
DefAns.List = 'Green';
DefAns.ListRank = 5;
DefAns.Result = {'val1',1,10,false;'val2',0,16,true};
DefAns.Vector = [1 3 6];
DefAns.VectorStrings = {'val1';'val2'};
DefAns.ListVector = [1 3 6];
DefAns.ListVectorStrings = {'val1';'val2'};
DefAns.Site = {'www.auction1.com','www.bidme.com'};
DefAns.SiteRank = [1 2 4];
DefAns.Date = [2013 1 2 0 0 0];
DefAns.Color = [1 0 0];
DefAns.X = 123584;
DefAns.Y = 94531;
DefAns.Z = [];
DefAns.Button = 2;
DefAns.Image = uint8(256*ones(100,100,3));
DefAns.hidden = 'HiddenVariable';
clearvars d files

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
