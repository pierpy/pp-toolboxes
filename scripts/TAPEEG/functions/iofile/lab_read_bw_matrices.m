function matrixout = lab_read_bw_matrices(Filename,doavg)

if ~exist('doavg','var')
    doavg = true;
end

if ~exist('Filename','var') | isempty(Filename)
   [Filename,Filepath]=uigetfile('*.*','Select file');
   Matrix_File = fullfile(Filepath,Filename);
else
   Matrix_File = Filename;
end

% load matrix-file
matrixout=[];
fid=fopen(Matrix_File);
tline = fgetl(fid);
if length(tline) < 5 | ~strcmp(tline(1:5),'File:')
    % try simple text-file import (eg BrainWave export or TAPEEG txt-export)
    matrixdata = [];
    while ischar(tline)
        matrixdata = [matrixdata;str2num(tline)]; %#ok<ST2NM,AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);
    if size(matrixdata,1) == size(matrixdata,2)
        matrixout.matrix = matrixdata;
        [~,~,~,FilenameS] = lab_filename(Matrix_File);
        matrixout.name = cellstr(FilenameS);
    else
        matrixout = [];
    end
elseif strcmp(tline(1:5),'File:') % = BrainWave Matrices from batch processing
    nummatrix = 0;
    doepoch = false;
    while ischar(tline)
        if length(tline) >= 5 & strcmp(tline(1:5),'File:')
            nummatrix = nummatrix + 1;
            if length(tline) > 5
                [~,~,~,name] = lab_filename(strtrim(tline(6:end)));
            else
                name = ['Matrix_' num2str(nummatrix)];
            end
        end
        if length(tline) >= 6 & strcmp(tline(1:6),'Epoch:')
            if length(tline) > 6
                epoch = str2num(tline(7:end)); %#ok<ST2NM>
                if epoch > 1
                    doepoch = true;
                end
            else
                epoch = 1;
            end
        end
        if ~isempty(tline) & ~isnan(str2double(tline(1)))
            matrixdata = sscanf(strrep(tline,',','.'),'%f')';
            for i = 2:size(matrixdata,2)
                tline = fgetl(fid);
                tmp = sscanf(strrep(tline,',','.'),'%f')';
                if size(tmp,2) == size(matrixdata,2)
                    matrixdata = [matrixdata;tmp]; %#ok<AGROW>
                else
                    disp('    Abort, error reading Brainwave Matrix-file')
                    matrixout = [];
                    return
                end
            end
            matrixout.matrix(:,:,nummatrix) = matrixdata;
            name1{nummatrix} = name; %#ok<AGROW>
            name2{nummatrix} = [name '_Epoch' num2str(epoch)]; %#ok<AGROW>
        end
        tline = fgetl(fid);  
    end
    fclose(fid);
    
    if exist('name2','var') & doepoch == true
        matrixout.name = name2;
    elseif exist('name1','var')
        matrixout.name = name1;
    else
        matrixout = [];
        return
    end
    
    if doavg == true
        prompt = matrixout.name{1};
        tmp = union(strfind(prompt,'_'),strfind(prompt,' '));
        if ~isempty(tmp)
            if tmp(1) > 1
                subjectname = prompt(1:tmp(1)-1);
            else
                subjectname = [];
            end
            for i = 1:length(tmp)
                tmp2 = num2str(i-1);
                numi = ' ';
                for j = 1:length(tmp2)
                    numi = [numi '^' tmp2(j)]; %#ok<AGROW>
                end
                if i ~= length(tmp)
                    subjectname = [subjectname numi ' ' prompt(tmp(i)+1:tmp(i+1)-1)]; %#ok<AGROW>
                elseif tmp(i) == length(prompt)
                    subjectname = [subjectname numi]; %#ok<AGROW>
                else
                    subjectname = [subjectname numi ' ' prompt(tmp(i)+1:end)]; %#ok<AGROW>
                end
            end
        else
            subjectname = prompt;
        end
        clearvars prompt tmp tmp2 numi
        
        settings.subjectname = 0;
        settings.AVGsubject = true;
        
        Prompt = cell(0,2);
        Formats = {};
        
        Prompt(end+1,:) = {'Average matrices per subject','AVGsubject'};
        Formats(end+1,1).type = 'check';
        Formats(end,1).size = [-1 -1];
        
        Prompt(end+1,:) = {'',''};
        Formats(end+1,1).type = 'text';
        
        Prompt(end+1,:) = {'Number of underscores in subject name',''};
        Formats(end+1,1).type = 'text';
        
        Prompt(end+1,:) = {['(' subjectname ')'],'subjectname'};
        Formats(end+1,1).type = 'edit';
        Formats(end,1).format = 'integer';
        Formats(end,1).size = 40;
        Formats(end,1).limits = [-1 20];
        
        [settings,Cancelled] = inputsdlg(Prompt,'Import BrainWave',Formats,settings);
        if Cancelled == 1
            return
        end
        
        if settings.AVGsubject == true
            subjectsrecall = '';
            subjectsnr = 0;
            for j = 1:nummatrix
                tmp = strfind(matrixout.name{j},'_');
                if size(tmp,2) > settings.subjectname
                    subjectstmp = matrixout.name{j}(1:tmp(settings.subjectname + 1)-1);
                elseif ~isempty(tmp)
                    subjectstmp = matrixout.name{j}(1:tmp(end)-1);
                else
                    subjectstmp = matrixout.name{j};
                end
                if strcmp(subjectstmp,subjectsrecall)
                    subjectstrials{subjectsnr} = [subjectstrials{subjectsnr} j]; %#ok<AGROW>
                else
                    subjectsnr = subjectsnr + 1;
                    subjectstrials{subjectsnr} = j; %#ok<AGROW>
                    subjects{subjectsnr} = subjectstmp; %#ok<AGROW>
                    subjectsrecall = subjectstmp;
                end
            end
            for i = 1:size(subjectstrials,2)
                matrixtmp(:,:,i) = mean(matrixout.matrix(:,:,subjectstrials{i}),3); %#ok<AGROW>
            end
            matrixout.matrix = matrixtmp;
            matrixout.name = subjects;
        end
    end
else
    matrixout = [];
end
