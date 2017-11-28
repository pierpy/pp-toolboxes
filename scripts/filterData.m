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
    nome_file=files(k).name
    cd(nome_file);
    conditions = dir;
    conditions(1:2)=[];
    conditions_path = strcat(files_path, '\', nome_file);
     for i = 1:1 %size(conditions, 1)
        condition_name = conditions(i).name
        cd(condition_name);
        periods = dir; 
        periods(1:2)=[];
        period_path = strcat(conditions_path, '\', condition_name);
            for j = 1:1% size(periods, 1)
                period_name = periods(j).name
                cd(period_name);
                txtFiles=dir('*.txt');
                txt_files_path = strcat(period_path, '\', period_name);
                  for n = 1: size(txtFiles, 1)
                        segment = load(strcat(txt_files_path,'\',txtFiles(n).name));
%                          segment=segment';
%                         [smoothdata] = eegfilt(segment,256,1,40);
% %                       [smoothdata] =ft_preproc_bandpassfilter(segment, 256, [1 40]);
                      saveeph(strcat(txt_files_path,'\',txtFiles(n).name(1:end-7),'_interp.ep'),segment);
                        clear segment smoothdata  ;
                  end
                  cd(period_path);
            end
            cd(conditions_path);
     end
     cd(files_path);
end


