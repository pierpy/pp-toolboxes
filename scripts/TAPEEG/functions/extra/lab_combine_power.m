% Script to import power results in different frequency bands, whereas
% results of every frequency band is stored in a different file. Output is
% a xls file, stored at the same location as the last input file
%
% written by F. Hatz 2012

function lab_combine_power

skipprocessing = 0;

data = [];
dataM = [];
stopreading = 0;
readNr = 0;
while stopreading == 0
    if skipprocessing == 0
        [filename,filepath]=uigetfile('*.*',['Select xls-file Nr ' num2str(readNr+1) '(absolute mean values)']);
        if ~isempty(filename) & filename ~= 0
            readNr = readNr + 1;
            cd(filepath);
            if ~exist('fileoutput','var')
                [~,~,~,filenameS] = lab_filename(filename);
                fileoutput = fullfile(filepath,filenameS);
            end
            if ispc
                [~,~,datatmp] = xlsread(fullfile(filepath,filename));
            else
                [~,~,datatmp] = xlsread(fullfile(filepath,filename),1,'','basic');
            end
            % find number of clustervars
            if ~exist('clustervars','var')
                for i = 2:size(datatmp,1)
                    tmp = strfind(datatmp{i,1},'_');
                    tmp2{i-1,1} = datatmp{i,1}(1:tmp-1);
                end
                [~,m,~]=unique(tmp2,'last');
                m = sort(m);
                n = m - m(1);
                if m(1:end-1) == n(2:end)
                    clustervars = m(1);
                end
                clearvars tmp tmp2 m n i
                clustervars = inputdlg('Number of variables in cluster','Nr variables',[1 50],{num2str(clustervars)});
                if ~isempty(clustervars)
                    clustervars = str2num(clustervars{1,1});
                else
                    skipprocessing = 1;
                end
            end
        else
            stopreading = 1;
        end
    end
    if skipprocessing == 0 & stopreading == 0
        if (clustervars * readNr)+1 == size(datatmp,1)
            stopreading = 1;
        end
        if isempty(data)
            header = datatmp(1,:);
        end
        data = cat(3,data,datatmp((readNr-1)*clustervars+2:readNr*clustervars+1,:));
        clearvars datatmp
        if exist(fullfile(filepath,[filename(1:end-6) 'dian.xls']),'file')
            if ispc
                [~,~,datatmp] = xlsread(fullfile(filepath,[filename(1:end-6) 'dian.xls']));
            else
                [~,~,datatmp] = xlsread(fullfile(filepath,[filename(1:end-6) 'dian.xls']),1,'','basic');
            end
            if isempty(dataM)
                headerM = datatmp(1,:);
            end
            dataM = cat(3,dataM,datatmp((readNr-1)*clustervars+2:readNr*clustervars+1,:));
            clearvars datatmp
        end
    end
end
if skipprocessing == 0
    if ~isempty(data) & exist('header','var')
        dataLog = data;
        dataR = data;
        dataRLog = data;
        datatmp = cell2mat(data(:,2:end,:));
        dataLog(:,2:end,:) = num2cell(log(datatmp));
        totpower = sum(datatmp,3);
        datatmp = datatmp ./ repmat(totpower,[1 1 size(datatmp,3)]);
        dataR(:,2:end,:) = num2cell(datatmp);
        dataRLog(:,2:end,:) = num2cell(lab_logit(datatmp));
        data = permute(data,[1 3 2]);
        data = reshape(data,size(data,1)*size(data,2),size(data,3));
        data = cat(1,header,data);
        if size(data,2) > 255
            fileout = [fileoutput '_Comb.xlsx'];
        else
            fileout = [fileoutput '_Comb.xls'];
        end
        lab_write_xls(fileout,data);
        
        dataLog = permute(dataLog,[1 3 2]);
        dataLog = reshape(dataLog,size(dataLog,1)*size(dataLog,2),size(dataLog,3));
        dataLog = cat(1,header,dataLog);
        if size(dataLog,2) > 255
            fileout = [fileoutput '_Comb_Log.xlsx'];
        else
            fileout = [fileoutput '_Comb_Log.xls'];
        end
        lab_write_xls(fileout,dataLog);
        
        dataR = permute(dataR,[1 3 2]);
        dataR = reshape(dataR,size(dataR,1)*size(dataR,2),size(dataR,3));
        dataR = cat(1,header,dataR);
        if size(dataR,2) > 255
            fileout = [fileoutput '_Comb_Realtive.xlsx'];
        else
            fileout = [fileoutput '_Comb_Realtive.xls'];
        end
        lab_write_xls(fileout,dataR);
        
        dataRLog = permute(dataRLog,[1 3 2]);
        dataRLog = reshape(dataRLog,size(dataRLog,1)*size(dataRLog,2),size(dataRLog,3));
        dataRLog = cat(1,header,dataRLog);
        if size(dataRLog,2) > 255
            fileout = [fileoutput '_Comb_RealtiveLogit.xlsx'];
        else
            fileout = [fileoutput '_Comb_RealtiveLogit.xls'];
        end
        lab_write_xls(fileout,dataRLog);
    end
    if ~isempty(dataM) & exist('headerM','var')
        dataMLog = dataM;
        dataMR = dataM;
        dataMRLog = dataM;
        datatmp = cell2mat(dataM(:,2:end,:));
        dataMLog(:,2:end,:) = num2cell(log(datatmp));
        totpower = sum(datatmp,3);
        datatmp = datatmp ./ repmat(totpower,[1 1 size(datatmp,3)]);
        dataMR(:,2:end,:) = num2cell(datatmp);
        dataMRLog(:,2:end,:) = num2cell(lab_logit(datatmp));
        dataM = permute(dataM,[1 3 2]);
        dataM = reshape(dataM,size(dataM,1)*size(dataM,2),size(dataM,3));
        dataM = cat(1,headerM,dataM);
        if size(dataM,2) > 255
            fileout = [fileoutput(1:end-2) 'dian_Comb.xlsx'];
        else
            fileout = [fileoutput(1:end-2) 'dian_Comb.xls'];
        end
        lab_write_xls(fileout,dataM);
        dataMLog = permute(dataMLog,[1 3 2]);
        dataMLog = reshape(dataMLog,size(dataMLog,1)*size(dataMLog,2),size(dataMLog,3));
        dataMLog = cat(1,headerM,dataMLog);
        if size(dataMLog,2) > 255
            fileout = [fileoutput(1:end-2) 'dian_Comb_Log.xlsx'];
        else
            fileout = [fileoutput(1:end-2) 'dian_Comb_Log.xls'];
        end
        lab_write_xls(fileout,dataMLog);
        dataMR = permute(dataMR,[1 3 2]);
        dataMR = reshape(dataMR,size(dataMR,1)*size(dataMR,2),size(dataMR,3));
        dataMR = cat(1,headerM,dataMR);
        if size(dataMR,2) > 255
            fileout = [fileoutput(1:end-2) 'dian_Comb_Relative.xlsx'];
        else
            fileout = [fileoutput(1:end-2) 'dian_Comb_Relative.xls'];
        end
        lab_write_xls(fileout,dataMR);
        dataMRLog = permute(dataMRLog,[1 3 2]);
        dataMRLog = reshape(dataMRLog,size(dataMRLog,1)*size(dataMRLog,2),size(dataMRLog,3));
        dataMRLog = cat(1,headerM,dataMRLog);
        if size(dataMRLog,2) > 255
            fileout = [fileoutput(1:end-2) 'dian_Comb_RelativeLogit.xlsx'];
        else
            fileout = [fileoutput(1:end-2) 'dian_Comb_RelativeLogit.xls'];
        end
        lab_write_xls(fileout,dataMRLog);
    end
end
