function [subjfullpath, correctname] = itab_createfolderstr( fullfilepath, analysis_folder_path, sbjcode, condcode, sesscode, reccode )
%ITAB_CREATEFOLDERSTR Summary of this function goes here
%   Detailed explanation goes here
    % alldatapath, sbjcode, condcode, sesscode, reccode
    if nargin == 1
        analysis_folder_path = uigetdir([],'Provide the path of the folder in wich do analysis (.edf files)');
        prompt = {'Enter subject code (6 letters):', 'Enter condition code (2 letters):', 'Enter session code (2 numbers e.g. 01):' , 'Enter record code (2 numbers e.g. 01):'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'xxxxxx', 'CC', '01', '01'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        correctname = strcat(answer{1}, answer{2}, '_', answer{3}, answer{4});
    elseif nargin == 2 || nargin ~= 5
        prompt = {'Enter subject code (6 letters):', 'Enter condition code (2 letters):', 'Enter session code (2 numbers e.g. 01):' , 'Enter record code (2 numbers e.g. 01):'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'xxxxxx', 'CC', '01', '01'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        correctname = strcat(answer{1}, answer{2}, '_', answer{3}, answer{4});
    elseif nargin == 5    
        correctname = strcat(sbjcode, condcode, '_', sesscode, reccode);
    else
        error('filename and analysis folder path are required');
    end
    [~, ~, ext] = fileparts(fullfilepath);
    subjfold = correctname(1:end-2);
    subjfullpath = strcat(analysis_folder_path, '\', subjfold);
    mkdir(subjfullpath)
    newfilepath = fullfile(analysis_folder_path, subjfold, strcat(correctname, ext));
    copyfile(fullfilepath, newfilepath)
end

