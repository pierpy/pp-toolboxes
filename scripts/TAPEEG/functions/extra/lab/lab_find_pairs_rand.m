% Helper-file for lab_find_pairs
%
% written by F. Hatz Vumc 2013

function [Pdata,Cdata,Repeats,ResultVMean,ResultVStd] = lab_find_pairs_rand(Pdata,Cdata,varsT,mode)

skipprocessing = 0;

tmp = inputdlg({'Stop after ... repeats','Max sequences'},'',[1 30;1 30],{'4','100000'});
if ~isempty(tmp)
    Repeats = str2num(tmp{1,1}); %#ok<ST2NM>
    maxcounter = str2num(tmp{2,1}); %#ok<ST2NM>
else
    skipprocessing = 1;
end
clearvars tmp
pause(0.2);

if skipprocessing == 0
    Pdat = cell2mat(Pdata(:,2:end));
    Cdat = cell2mat(Cdata(:,2:end));
    for i = 1:size(varsT,2)
        Pdat(:,i) = Pdat(:,i)/varsT(i);
        Cdat(:,i) = Cdat(:,i)/varsT(i);
    end
end

ResultVMean = [];
counter = 0;
counterall = 0;
%Diffs = sum(abs(repmat(permute(Pdat,[1 3 2]),[1 size(Cdat,1) 1]) - repmat(permute(Cdat,[1 3 2])',[size(Pdat,1) 1 1])),3);

progressbar;
progressbar('Find pairs by random approach');

while max(counter) < Repeats & counterall < maxcounter
    tmp = randperm(size(Pdat,1));
    Pdattmp = Pdat(tmp,:);
    RtmpP = Pdata(tmp,:);
    clearvars tmp
    Cdattmp = Cdat;
    Cdatatmp = Cdata;
    Ccounter = (1:size(Cdattmp,1))';
    skipturn = 0;
    for j = 1:size(Pdat,1)
        if skipturn == 0
            diffs = abs(Cdattmp - repmat(Pdattmp(j,:),size(Cdattmp,1),1));
            if strcmp(mode,'error-to-fit')
                diffs = sum(diffs,2);
                tmp = find(diffs == min(diffs),1,'first');
            else
                tmp = 1:size(diffs,1);
                for i = 1:size(Pdat,2)
                    tmp2 = find(diffs(tmp,1) < varsT(i));
                    if isempty(tmp2) & i == 1
                        skipturn = 1;
                        disp('   Abort: matching not possible')
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
                clearvars varN tmp2
                diffs = sum(diffs,2);
            end
            RtmpC(j,:) = Cdatatmp(tmp,:); %#ok<AGROW>
            RtmpV(j,1) = diffs(tmp); %#ok<AGROW>
            RtmpCounter(j,1) = Ccounter(tmp); %#ok<AGROW>
            tmp = setdiff(1:size(Cdattmp,1),tmp);
            Cdatatmp = Cdatatmp(tmp,:);
            Cdattmp = Cdattmp(tmp,:);
            Ccounter = Ccounter(tmp,1);
            clearvars tmp diffs
        end
    end
    if skipturn == 0
        ResultMean = mean(RtmpV);
        ResultStd = std(RtmpV);
        datatmp = sortrows([RtmpP RtmpC num2cell(RtmpCounter)],1);
        RtmpP = datatmp(:,1:size(RtmpP,2));
        RtmpC = datatmp(:,size(RtmpP,2)+1:end-1);
        RtmpCounter = cell2mat(datatmp(:,end));
        clearvars datatmp
        if isempty(ResultVMean) | (ResultVMean+ResultVStd) > (ResultMean+ResultStd)
            ResultVMean = ResultMean;
            ResultVStd = ResultStd;
            clearvars ResultB ResultC
            ResultC{1,1} = RtmpC;
            ResultP = RtmpP;
            ResultB(:,1) = RtmpCounter;
            counter = 1;
            disp([num2str(ResultVMean) '+-' num2str(ResultVStd) ': first'])
        %elseif ResultVMean == ResultMean & ResultVStd > ResultStd
        %    ResultVMean = ResultMean;
        %    ResultVStd = ResultStd;
        %    clearvars ResultB ResultC
        %    ResultC{1,1} = RtmpC;
        %    ResultP = RtmpP;
        %    ResultB(:,1) = RtmpCounter;
        %    counter = 1;
        %    disp([num2str(ResultVMean) '+-' num2str(ResultVStd) ': first'])
        elseif (ResultVMean+ResultVStd) == (ResultMean+ResultStd)
            foundequal = 0;
            for j = 1:size(ResultB,2)
                if isequal(ResultB(:,j),RtmpCounter)
                    counter(j) = counter(j) + 1; %#ok<AGROW>
                    foundequal = 1;
                end
            end
            if foundequal == 0
                ResultC{1,end+1} = RtmpC; %#ok<AGROW>
                ResultB(:,end+1) = RtmpCounter; %#ok<AGROW>
                counter(end+1) = 1; %#ok<AGROW>
            end
            disp([num2str(ResultVMean) '+-' num2str(ResultVStd) ': ' num2str(max(counter)) ' repeats (variations: ' num2str(size(ResultB,2)) ')'])
        end
        clearvars RtmpC RtmpP RtmpB RtmpV Pdattmp Cdattmp Cdatatmp Resulttmp
    end
    counterall = counterall + 1;
    if floor(counterall/1000) == counterall/1000
        progressbar(counterall/maxcounter);
    end
end
clearvars i
Cdata = ResultC{1,counter == max(counter)};
Pdata = ResultP;
Repeats = max(counter);
progressbar(1);

end