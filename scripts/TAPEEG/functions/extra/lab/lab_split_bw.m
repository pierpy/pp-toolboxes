function lab_split_bw

[filename,filepath] = uigetfile('*.*','Select BrainWave file');

[xlsfilename,xlsfilepath] = uigetfile('*.*','Select file');
[results,subjectstmp] = xlsread(fullfile(xlsfilepath,xlsfilename));
if size(subjectstmp,1) < size(subjectstmp,2)
    subjectstmp = subjectstmp';
    results = results';
end
strlist = subjectstmp(1,2:end);
selection = listdlg('PromptString','Select result','SelectionMode','single','ListString',strlist);
results = results(:,selection);
resultsnr = unique(results);
for i = 1:length(resultsnr)
    subjects{i,1} = subjectstmp(find(results == resultsnr(i))+1,1);
end

for i = 1:size(subjects,1)
    fid=fopen(fullfile(filepath,filename),'r');
    fidout=fopen(fullfile(filepath,[filename(1:end-4) '_' num2str(i) '.txt']),'w');
    tline = fgets(fid);
    fprintf(fidout,tline);
    tline = fgets(fid);
    while ~isnumeric(tline)
        if length(tline) > 10
            tmp = strfind(tline,'_');
            tmp = tline(1:tmp(1)-1);
            tmp2 = regexp(tmp,'\d');
            if length(tmp2) == length(tmp)
                tmp = ['P_' tmp];
            end
            if max(strcmp(subjects{i,1},tmp))
                fprintf(fidout,tline);
                nobreak = 0;
            else
                nobreak = 1;
            end
        else
            if nobreak == 0
                fprintf(fidout,tline);
            end
        end
        tline = fgets(fid);
    end
    fclose(fid);
    fclose(fidout);
end