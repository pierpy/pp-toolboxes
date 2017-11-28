    function [res,Names,msg] = SlurpData(path,conds,directory, tp,pb)
% This function helps you to read well organized data from a directory into matlab
% PATH is a wildcart to identify the files of one condition across all
% subjects in the directory
% CONDS is a cell array containing the tags for all conditions you want to
% select. The condition used in the path must be contained in this cell
% array
%
% The tool will go across all subjects found with the wildcart and read all
% conditions in the sequence defined by CONDS
%
% The output is a 4D Time by Condition by Subject by Electrode matrix
% The optional additional output is a cell array with the filenames read
%
% Have fun
%
% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

if (nargin < 3)
    directory = pwd();
end

if (nargin < 4)
    tp = 0;
end


if (nargin < 5)
    pb = 0;
end

d = dir(fullfile(directory,path));

ds = struct2cell(d);
NameList = ds(1,:);
NameListSorted = lower(sortrows(NameList'));

conds = lower(conds);


if numel(d) == 0
    if nargout > 2
        msg = 'No file found';
        res = [];
        Names = [];
        return
    else
        error('No file found');
    end
end

idx = -1;
for j = 1:numel(conds)
    if ~isempty(strfind(lower(d(1).name),conds{j}))
        idx = j;
    end
end

if idx == -1
    if nargout > 2
        msg = 'No condition tag found in the choosen path';
        res = [];
        Names = [];
        return;
    else
        error('No condition tag found in the choosen path');
    end
end

if pb == 1
    h = waitbar(0);
    set(h,'Name','Importing ASCII data');
end

Names = cell(numel(d),numel(conds));
cnt = 1;
for i = 1:numel(d)
    if numel(strfind(NameListSorted{i},conds{idx})) ~= 1
        msg = ['No or several condition tags found in ' NameListSorted{i}];
        res = [];
        Names = [];
        if (pb == 1)
            close(h);
        else
            error(msg);
        end
        return;
    end
    for j = 1:numel(conds)
        
        fn = strrep(NameListSorted{i},conds{idx},conds{j});
       
        if (exist(fullfile(directory,fn),'file') == 0)
            res = [];
            Names = [];
            msg = ['Dataset was not complete, import aborted, file ' fn ' (and maybe others) missing'];

            if (pb == 0)
                error([fn ' does not exist']);
            else
                close(h);
            end
            return;
        end
       
        if (pb == 0)
            disp(fn);
        else
            waitbar(cnt / numel(d) / numel(conds),h);
            set(get(get(h,'Children'),'Title'),'String',sprintf('Importing %s',fn),'Interpreter','none');
        end
        
        
        cnt = cnt +1;
            
        if ((i == 1) && (j == 1)) && (tp == 1)
            res(i,j,:,:) = load(fullfile(directory,fn));
            Names{i,j} = fn;
        elseif ((i == 1) && (j == 1)) && (tp == 0)
            res(i,j,:,:) = load(fullfile(directory,fn))';
            Names{i,j} = fn;
            
        else
            if (tp == 1)
                    tmp = load(fullfile(directory,fn));
            else
                    tmp = load(fullfile(directory,fn))';
            end

            if (size(tmp,1) ~= size(res,3) || size(tmp,2) ~= size(res,4))
                if nargout > 2
                    msg = ['Sample size not constant. File causing the problem is ' fn];
                    res = [];
                    Names = [];
                    return;
                else
                    error('Sample size not constant');
                end
            end

            res(i,j,:,:) = tmp;
            Names{i,j} = fn;
        end
    end
end

msg = sprintf('%i observations with %i condition loaded (%i frames x %i channels)',size(res,1),size(res,2),size(res,4),size(res,3));

if nargout < 3
    disp(msg);
end

if pb == 1
    close(h);
end