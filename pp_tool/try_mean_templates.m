clear
close all
clc
path='/Users/pp/Documents/MATLAB/microstates/dati/mapsControlli';

cd(path);
files=dir('*.ep');
maps=template;


Nm = 4;

for k=1:size(files,1)
    segment = load(strcat(path,'/',files(k).name));
    maps=[maps;segment];   
end

[final_templates, stab, V]=mean_templates(maps, Nm);
plot(final_templates')