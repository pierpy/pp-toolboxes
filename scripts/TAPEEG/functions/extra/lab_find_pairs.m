% Script to find matching pairs in two different groups
%
% Result = lab_find_pairs(data)
%          data: filepath/filename to xls-file
%                cellmatrix with xls-content
%
% Input is a xls-file with 1 patient per row, first row must contain header
% Group information must be one of the variables, lower number patient
% group, at least same number of control patients
%
% Output is xls-file with paired subjects
%
% written by F. Hatz 2013

function Result = lab_find_pairs

disp('Find pairs in two groups (e.g. patients and healthy controls)')

Result = [];
[data,header,result,~,cfg] = lab_read_statistics([],1,0,0,1,1);
if isempty(data)
    return
end
result = result(:,1);
[result,Ridx] = sort(result);
data = data(Ridx,:);
header.subjects = header.subjects(Ridx);
[~,tmp] = unique(result,'last');
if length(tmp) < 2
    disp('Abort: invalid grouping variable')
    return
end
Pdata = data(1:tmp(1),:);
Psubj = header.subjects(1:tmp(1));
Cdata = data(tmp(1)+1:tmp(2),:);
Csubj = header.subjects(tmp(1)+1:tmp(2));
if size(Pdata,1) > size(Cdata,1)
    tmp = Pdata;
    Pdata = Cdata;
    Cdata = tmp;
    tmp = Psubj;
    Psubj = Csubj;
    Csubj = tmp;
    clearvars tmp
end

FORMAT = [];
settings = [];
for i = 1:size(data,2)
    settings(i).Name = header.vars{i}; %#ok<AGROW>
    if min(data(:,i)) == 0 & max(data(:,i)) == 1
        settings(i).Threshold = 1; %#ok<AGROW>
    else
        settings(i).Threshold = abs((max(data(:,i)) - min(data(:,i))) / 5); %#ok<AGROW>
    end
    settings(i).MinValue = num2str(min(data(:,i))); %#ok<AGROW>
    settings(i).MaxValue = num2str(max(data(:,i))); %#ok<AGROW>
    FORMAT{strcmp(fieldnames(settings),'Name')} = 'text'; %#ok<AGROW>
    FORMAT{strcmp(fieldnames(settings),'MinValue')} = 'text'; %#ok<AGROW>
    FORMAT{strcmp(fieldnames(settings),'MaxValue')} = 'text'; %#ok<AGROW>
    tmp = [];
    tmp.type = 'edit';
    tmp.format = 'float';
    tmp.size = 50;
    tmp.limits = [-inf inf];
    FORMAT{strcmp(fieldnames(settings),'Threshold')} = tmp; %#ok<AGROW>
    clearvars tmp
end
[settings,skipprocessing] = inputsdlg(settings,'Thresholds',FORMAT);
if skipprocessing == 1
    return
end
varsT = zeros(1,length(settings));
for i = 1:length(settings)
    varsT(i) = settings(i).Threshold;
end
clearvars settings

% set mode
settings.mode = 'sequential';
settings.dorandom = false;
Prompt = cell(0,2);
Formats = [];

Prompt(end+1,:) = {'Mode','mode'};
Formats(end+1,1).type = 'list';
Formats(end,1).style = 'popupmenu';
Formats(end,1).format = 'input'; % Answer will give value shown in items, disable to get integer
Formats(end,1).items = {'sequential','error-to-fit'};

Formats(end+1,1).type = 'none';

Prompt(end+1,:) = {'Randomized sequence search','dorandom'};
Formats(end+1,1).type = 'check';

[settings,skipprocessing] = inputsdlg(Prompt,'Mode',Formats,settings);
if skipprocessing == 1
    return
end

if settings.dorandom == false
    % sort patients (1.upper extreme, 1. lower extreme, 2.upper extreme...)
    vars = 1:size(Pdata,2);
    Pdata2 = [Psubj num2cell(Pdata)];
    Pdata2 = sortrows(Pdata2,-vars);
    Pdata = cell2mat(Pdata2(:,2:end));
    Psubj = Pdata2(:,1);
    tmp = 1:size(Pdata,1);
    tmp2(1:2:size(tmp,2)) = tmp(1:round(size(tmp,2)/2));
    if floor(size(tmp,2)/2) ~= size(tmp,2)/2
        tmp2(size(tmp,2)-1:-2:1) = tmp(end-floor(size(tmp,2)/2)+1:end);
    else
        tmp2(size(tmp,2):-2:1) = tmp(end-floor(size(tmp,2)/2)+1:end);
    end
    Pdata = Pdata(tmp2,:);
    Psubj = Psubj(tmp2);
    clearvars tmp tmp2
end
% extract variable data
for i = 1:size(varsT,2)
    Pdata(:,i) = Pdata(:,i)/varsT(i);
    Cdata(:,i) = Cdata(:,i)/varsT(i);
end

if settings.dorandom == false
    Cdattmp = Cdata;
    Csubjtmp = Csubj;
    for j = 1:size(Pdata,1)
        diffs = abs(Cdattmp - repmat(Pdata(j,:),size(Cdattmp,1),1));
        if strcmp(settings.mode,'error-to-fit')
            diffs = sum(diffs,2);
            tmp = find(diffs == min(diffs),1,'first');
        else
            tmp = 1:size(diffs,1);
            for i = 1:size(vars,2)
                tmp2 = find(diffs(tmp,1) < 1);
                if isempty(tmp2) & i == 1
                    disp('   Abort: matching not possible')
                    return
                elseif ~isempty(tmp2)
                    tmp = tmp(tmp2);
                end
            end
            varN = 1;
            while length(tmp) > 1
                tmp2 = diffs(tmp,varN) == min(diffs(tmp,varN));
                tmp = tmp(tmp2);
                varN = varN + 1;
                if varN > size(diffs,2)
                    tmp = tmp(1);
                end
            end
        end
        Stmp(j,:) = Csubjtmp(tmp,:); %#ok<AGROW>
        diffs = sum(diffs,2);
        RtmpV(j,1) = diffs(tmp); %#ok<AGROW>
        tmp = setdiff(1:size(Cdattmp,1),tmp);
        Cdattmp = Cdattmp(tmp,:);
        Csubjtmp = Csubjtmp(tmp,:);
        clearvars tmp tmp2 diffs varN
    end
    Csubj = Stmp;
    Cdata = Cdattmp;
    ResultVMean = mean(RtmpV);
    ResultVStd = std(RtmpV);
else
    Pdata2 = [Psubj num2cell(Pdata)];
    Cdata2 = [Csubj num2cell(Cdata)];
    [Pdata2,Cdata2,Repeats,ResultVMean,ResultVStd] = lab_find_pairs_rand(Pdata2,Cdata2,varsT,settings.mode);
    if isempty(ResultVMean)
        return
    end
    Pdata = cell2mat(Pdata2(:,2:end));
    Psubj = Pdata2(:,1);
    Cdata = cell2mat(Cdata2(:,2:end));
    Csubj = Cdata2(:,1);
end

% Collect and write result
datatmp = sortrows([Pdata Cdata],1);
Pdata = datatmp(:,1:size(Pdata,2));
Cdata = datatmp(:,size(Pdata,2)+1:end);
clearvars datatmp
Pdata(:,end+1) = 0;
Cdata(:,end+1) = 1;
Bdata = [];
Bsubj = {};
Bdata(1:2:size(Pdata,1)*2,:) = Pdata;
Bdata(2:2:size(Pdata,1)*2+1,:) = Cdata;
Bsubj(1:2:size(Pdata,1)*2) = Psubj;
Bsubj(2:2:size(Pdata,1)*2+1) = Csubj;
xlsout = cat(2,{''},Bsubj(:)');
xlsout = cat(1,xlsout,cat(2,cat(1,header.vars(:),{'Group'}),num2cell(Bdata')));
[~,filepath,~,filenameS] = lab_filename(cfg.filename);
if size(Bdata,2) > 255
    fileout = fullfile(filepath,[filenameS '_paired.xlsx']);
else
    fileout = fullfile(filepath,[filenameS '_paired.xls']);
end
lab_write_xls(fileout,xlsout);

Result.all = Bdata(:,1:end-1);
Result.patients = Pdata(:,1:end-1);
Result.controls = Cdata(:,1:end-1);
Result.Pmean = mean(Pdata(:,1:end-1));
Result.Cmean = mean(Cdata(:,1:end-1));
Result.Pstd = std(Pdata(:,1:end-1));
Result.Cstd = std(Cdata(:,1:end-1));

% write verbose
fid=fopen(fullfile(filepath,[filenameS '_paired.vrb']),'w');
fprintf(fid,'Find pairs in 2 groups\n');
fprintf(fid,datestr(now,0));
fprintf(fid,'\n');
fprintf(fid,['Input File: ' filenameS '\n']);
fprintf(fid,['Ouput File: ' filenameS '_paired\n']);
fprintf(fid,'\n');
fprintf(fid,['Match criteria: ' settings.mode '\n']);
fprintf(fid,'\n');
if settings.dorandom == true
    fprintf(fid,['Randomized sequence search: Yes\n']);
    fprintf(fid,['   number of repeats: ' num2str(Repeats) ' (Error-to-fit: ' num2str(ResultVMean) '+-' num2str(ResultVStd) ')\n']);
else
    fprintf(fid,['Randomized sequence search: No\n']);
    fprintf(fid,['Error-to-fit: ' num2str(ResultVMean) '+-' num2str(ResultVStd) ')\n']);
end
fprintf(fid,'\n');
fprintf(fid,['Variables:     ' sprintf('%s|',header.vars{:}) '\n']);
fprintf(fid,['Bad threshold: ' num2str(varsT) '\n']);
fprintf(fid,'\n');
fprintf(fid,['Mean values (patients): ' num2str(Result.Pmean) ' (std:' num2str(Result.Pstd) ')\n']);
fprintf(fid,['Mean values (controls): ' num2str(Result.Cmean) ' (std:' num2str(Result.Cstd) ')\n']);
fclose(fid);

end
