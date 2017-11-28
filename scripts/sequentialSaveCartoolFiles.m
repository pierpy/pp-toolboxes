clear all
close all
clc

files_path='C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\_Microstate-Rest';

% prompt = 'cartella:';
% dlg_title = 'Input';
% str = input(prompt,'s');



% files_path=strcat(files_path,'\',str);
cd(files_path);
files=dir;
files(1:2)=[];

for k=1:1%size(files,1)
    nome_file=files(k).name;
    cd(nome_file);
    conditions = dir;
    conditions(1:2)=[];
    conditions_path = strcat(files_path, '\', nome_file);
     for i = 1: size(conditions, 1)
        condition_name = conditions(i).name;
        cd(condition_name);
        periods = dir; 
        periods(1:2)=[];
        period_path = strcat(conditions_path, '\', condition_name);
            for j = 1: size(periods, 1)
                period_name = periods(j).name;
                cd(period_name);
                txtFiles = dir;
                txtFiles(1:2)=[];
                txt_files_path = strcat(period_path, '\', period_name);
                  for n = 1: size(txtFiles, 1)
                        thedata=load(txtFiles(n).name);
                        filename = txtFiles(n).name;
                        savefilename=strcat(txt_files_path,'\',filename(1:end-7), 'ep');
%                         savefilename=strcat(savefilename,'ep');
                        saveeph(savefilename,thedata)
                        clear thedata filename;
                  end
                  cd(period_path);
            end
            cd(conditions_path);
     end
     cd(files_path);
end


