% demo script for inputsdlg_struct
%
% inputsdlg with struct as input allows to edit structure variables with size > 1

clear; close all;

% create structure variable with size > 1, DefAns = values, FORMATS = format of value if needed
for datanr = 1:5
    % no editable text / FORMAT = 'text'
    DefAns(datanr).Name = ['Data ' num2str(datanr)];
    FORMATS{strcmp(fieldnames(DefAns),'Name')} = 'text'; 
    
    % numeric vector / FORMAT = 'vector' 
    DefAns(datanr).Vector = [1 8 31];
    FORMATS{strcmp(fieldnames(DefAns),'Vector')} = 'vector'; 
    
    % numeric single value, no FORMAT needed
    DefAns(datanr).Value = -5; 
    
    % editable text, no FORMAT needed
    DefAns(datanr).String = 'string';
    
    % popupmenu with strings / FORMAT = 'cellarray with strings'
    DefAns(datanr).List = 'bluered';
    FORMATS{strcmp(fieldnames(DefAns),'List')} = {'bluered','autumn','bone','colorcube','cool','copper','gray','hot','hsv','jet'};
    
    % color, no FORMAT needed
    DefAns(datanr).Color = rand(1,3);
    
    % logical value, no FORMAT needed
    DefAns(datanr).AddPlot = false;
    
    % image, no Format needed
    DefAns(datanr).Image = uint8(256*ones(100,100,3)); % white image of size 100*100
end
clearvars datanr

Answer = inputsdlg(DefAns,'Inputsdlg_Struct',FORMATS);