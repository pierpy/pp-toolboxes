clear all
close all
clc

Nt=500;
Nm = 4;
Nch = 62;

V=zeros(Nt, Nch);
a = -2;
b = 2;
A = (b-a).*1*randn(Nm,Nt) + a;
T = (b-a).*1*randn(Nm,Nch) + a;

for i=1:50
eeg_or(i,:)=T(1,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=51:100
eeg_or(i,:)=T(2,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=101:150
eeg_or(i,:)=T(1,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=151:200
eeg_or(i,:)=T(4,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=201:250
eeg_or(i,:)=T(3,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=251:300
eeg_or(i,:)=T(2,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=301:350
eeg_or(i,:)=T(2,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=351:400
eeg_or(i,:)=T(2,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=401:450
eeg_or(i,:)=T(2,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end
for i=451:500
eeg_or(i,:)=T(1,:);
% eeg_or(i,:)=smooth(eeg_or(i,:));
end

eeg_or = eeg_or+0.01*rand(size(eeg_or));

%struttura con informazioni circa l'analisi
confstruct.substract_column_at_start = 1; % detrend si o no
confstruct.method_GFPeak = 'GFPL2';       % gfp calculation method
confstruct.use_gfp_peaks =0  ;             % utililla i picchi del gfp o tutto il file
confstruct.setall1 = 0;
confstruct.normalize = 0;
confstruct.tf = 500;
confstruct.fs = 500;
confstruct.nch = 62;
confstruct.minclusters = 4;
confstruct.maxclusters = 4;
confstruct.pmode = 1;
confstruct.algorithm = 1;
confstruct.similarity_measure = 'correlation';
confstruct.debug = 0;
confstruct.segment_min_len = 3;

[ eeg, gfp_peaks_indices, gfp_curve ] = modmaps_preprocess(eeg_or, confstruct);

[ NMicr ] = CW2( eeg ,gfp_peaks_indices, confstruct);

maps = NMicr.template{1, 1};
plot(gfp_curve)
[ output ] = compute_mstate_parameters( confstruct, eeg, maps );
