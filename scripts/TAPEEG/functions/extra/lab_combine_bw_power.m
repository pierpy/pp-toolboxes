% Script to import power results in different frequency bands, whereas
% results of every frequency band is stored in a different file. Output is
% a xls file, stored at the same location as the last input file
%
% written by F. Hatz 2012

function lab_combine_bw_power

skipprocessing = 0;

pnames = {'deltapower';'thetapower';'alpha1power';'alpha2power';'betapower';'gammapower'};
data = [];
for i = 1:6
    if skipprocessing == 0
        [filename,filepath]=uigetfile('*.*',['Select BrainWave results-file for ' pnames{i,1}]);
        if filename == 0
            return
        end
        cd(filepath);
        cfg.selection = [1 1+i];
        [datatmp,~,cfg] = lab_import_bw_results(fullfile(filepath,filename),cfg,'Single');
        if strcmp(datatmp{2,1}(1:10),'totalpower') & strfind(datatmp{cfg.clustervars+2,1},pnames{i})
            datatmp2 = cell2mat(datatmp(2:end,2:end));
            datatmp2 = datatmp2(1:cfg.clustervars,:) .* datatmp2(cfg.clustervars+1:end,:);
            if ~exist('result','var')
                result = datatmp;
            else
                result(end+1:end+cfg.clustervars,1) = datatmp(cfg.clustervars+2:end,1);
            end
            if size(datatmp,2) == size(result,2)
                data = cat(3,data,datatmp2);
                clearvars datatmp2
            else
                disp('wrong input file -- abort')
                skipprocessing = 1;
            end
        else
            disp('wrong input file -- abort')
            skipprocessing = 1;
        end
        clearvars datatmp
    end
end
totpower = sum(data,3);
data = data ./ repmat(totpower,[1 1 6]);
data = cat(3,totpower,data);
data = permute(data,[1 3 2]);
data = reshape(data,size(data,1)*size(data,2),size(data,3));
result(2:end,2:end) = num2cell(data);
[~,~,~,filename] = lab_filename(filename);
if size(result,2) > 255
    fileout = fullfile(filepath,['band' filename(6:end) '.xlsx']);
else
    fileout = fullfile(filepath,['band' filename(6:end) '.xls']);
end
lab_write_xls(fileout,result);

return
