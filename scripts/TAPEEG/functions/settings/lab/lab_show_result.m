% Helper-script for inputsdlg to show results

function result = lab_show_result(result)

if isempty(result)
    Fout = '';
    return
end
if isstruct(result)
    Fanswer = result;
    if length(Fanswer) > 1
        SizeFanswer = length(Fanswer);
        Fanswer = Fanswer(1);
    end
else
    Fanswer.value = result;
end
Ftmp = fieldnames(Fanswer);
for i = 1:length(Ftmp)
    if isempty(Fanswer.(Ftmp{i}));
        Fout{i,1} = [Ftmp{i} ':'];
    elseif size(Fanswer.(Ftmp{i}),1) > 1 & size(Fanswer.(Ftmp{i}),2) > 1;
        Fout{i,1} = [Ftmp{i} ': table' num2str(size(Fanswer.(Ftmp{i}),1)) 'x' num2str(size(Fanswer.(Ftmp{i}),2))];
    elseif isnumeric(Fanswer.(Ftmp{i})) | islogical(Fanswer.(Ftmp{i}));
        if size(Fanswer.(Ftmp{i})(1,:),2) > 3
            Fout{i,1} = [Ftmp{i} ': ' num2str(Fanswer.(Ftmp{i})(1,1:3)) '...'];
        else
            Fout{i,1} = [Ftmp{i} ': ' num2str(Fanswer.(Ftmp{i})(1,:))];
        end
    elseif ischar(Fanswer.(Ftmp{i}));
        if size(Fanswer.(Ftmp{i})(1,:),2) > 10
            Fout{i,1} = [Ftmp{i} ': ' Fanswer.(Ftmp{i})(1,1:10) '...'];
        else
            Fout{i,1} = [Ftmp{i} ': ' Fanswer.(Ftmp{i})(1,:)];
        end
    elseif iscell(Fanswer.(Ftmp{i}));
        if size(Fanswer.(Ftmp{i}),1) > 1 | size(Fanswer.(Ftmp{i}),2) > 1;
            Fout{i,1} = [Ftmp{i} ': table' num2str(size(Fanswer.(Ftmp{i}),1)) 'x' num2str(size(Fanswer.(Ftmp{i}),2))];
        elseif isnumeric(Fanswer.(Ftmp{i}){1});
            Fout{i,1} = [Ftmp{i} ': ' num2str(Fanswer.(Ftmp{i}){1})];
        elseif ischar(Fanswer.(Ftmp{i}){1})
            Fout{i,1} = [Ftmp{i} ': ' Fanswer.(Ftmp{i}){1}];
        else
            Fout{i,1} = [Ftmp{i} ': {' class(Fanswer.(Ftmp{i}){1}) '}'];
        end
    else
        Fout{i,1} = [Ftmp{i} ': {' class(Fanswer.(Ftmp{i})) '}'];
    end
end
if size(Fout,1) == 1
    Ftmp = strfind(Fout{1,1},':');
    if isempty(Ftmp) | Ftmp == length(Fout{1,1})
        Fout{1,1} = '';
    else
        Fout{1,1} = Fout{1,1}(Ftmp+2:end);
    end
end
if exist('SizeFanswer','var')
    Fout = cat(1,cellstr(['MULTIVARIABLE: ' num2str(SizeFanswer)]),Fout);
end
Fout = char(Fout);
clearvars Fanswer

uiwait(helpdlg(Fout,' '));