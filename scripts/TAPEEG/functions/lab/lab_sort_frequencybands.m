function [Filelist,Freqband] = lab_sort_frequencybands(Filelist)

Folders = cell(length(Filelist),1);
Filename = cell(length(Filelist),1);
for i = 1:length(Filelist)
    [Filename{i},Folders{i}] = lab_filename(Filelist{i});
end

dosort = true;
Folders2 = Folders;
Isort = [];
Freqs = zeros(length(Folders2),2);
while dosort == true
    Freqbands = zeros(length(Folders2),2);
    for i = 1:length(Folders2)
        if ~isempty(Folders2{i})
            Freqbands(i,:) = find_freqband(Folders2{i});
            if min(Freqbands(i,:)) > 0 & min(Freqs(i,:)) == 0
                Freqs(i,:) = Freqbands(i,:);
            end
        end
    end
    dosort = false;
    for i = 1:length(Folders2)
        [~,Folders2{i}] = lab_filename(Folders2{i});
        if ~isempty(Folders2{i})
            dosort = true;
        else
            Folders2{i} = ' ';
        end
    end
    [~,~,tmp1] = unique(Folders2,'stable');
    [tmp2,tmp3] = sortrows([tmp1(:) Freqbands],[1 2 3]);
    [~,~,tmp4] = unique(cellstr(num2str(tmp2)),'stable');
    tmp3(tmp3) = tmp4;
    Isort(:,end+1) = tmp3; %#ok<AGROW>
    clearvars tmp1 tmp2 tmp3 tmp4
end
Freqbands = zeros(length(Filename),2);
for i = 1:length(Filename)
    if ~isempty(Filename{i})
        [Freqbands(i,:),Filename{i}] = find_freqband(Filename{i});
        if min(Freqbands(i,:)) > 0 & min(Freqs(i,:)) == 0
            Freqs(i,:) = Freqbands(i,:);
        end
    end
end
[~,~,tmp1] = unique(Filename,'stable');
[tmp2,tmp3] = sortrows([tmp1(:) Freqbands],[1 2 3]);
[~,~,tmp4] = unique(cellstr(num2str(tmp2)),'stable');
tmp3(tmp3) = tmp4;
Isort(:,end+1) = tmp3;
clearvars tmp1 tmp2 tmp3 tmp4

[~,Isort2] = sortrows(Isort,size(Isort,2):-1:1);
Filelist = Filelist(Isort2);
Filelist = Filelist(:)';
Freqs = Freqs(Isort2,:);
Freqband = cell(1,size(Freqs,1));
for i = 1:size(Freqs,1)
    Freqband{1,i} = ['F' num2str(Freqs(i,1)) 'F' num2str(Freqs(i,2))];
end
    
end

function [freqband,input] = find_freqband(input)
    freqband = zeros(1,2);
    if isempty(input)
        return
    end
    tmp = strfind(input,'F');
    flag = false(1,length(tmp));
    for i = 1:length(tmp)
        if tmp(i)+1 <= length(input) & ~strcmp(input(tmp(i)+1),'i') & ~isnan(str2double(input(tmp(i)+1)))
            flag(1,i) = true;
        end
    end
    tmp = tmp(flag);
    if isempty(tmp)
        return
    end
    freqband(1,1) = sscanf(input(tmp(1)+1:end),'%d');
    if length(tmp) > 1
        freqband(1,2) = sscanf(input(tmp(2)+1:end),'%d');
    end
    input = regexprep(input,['F' num2str(freqband(1,1)) 'F' num2str(freqband(1,2))],'');
end