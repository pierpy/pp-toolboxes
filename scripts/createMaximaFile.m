function [] = createMaximaFile(subjectName, cartella)
% samplingF=256;
%subjectName='capger';
% path=pathString;
path=strcat('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiBetaMotorio\Dati_Microstati_MOTORIO\',subjectName,'\',cartella);
cd(path);
files=dir('*.interpolate.ep');
totalPeack=zeros(123,1);
for k=1:size(files,1)
    segment = load(strcat(path,'\',files(k).name));
    segment=segment';
    [peackData] = globalFieldPower(segment);
    numeroDiMassimi(k)=size(peackData,2);
    totalPeack=[totalPeack,peackData];   
end
totalPeack(:,1)=[];

saveeph(strcat(path,'\gfpMax_',cartella,subjectName,'.ep'),totalPeack');
% saveeph(strcat(path,'\gfpMax_',subjectName,'.eph'),totalPeack', samplingF);
