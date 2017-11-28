pathLeftAG_b='C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\_Microstate-Rest\Aleter\Left-AG\baseline';
cd(pathLeftAG_b);
elencoFiles=dir('*.Export_OK.INTERP.ep');
maps=load('C:\Users\Pierpaolo\Documents\MATLAB\MicrostatiStimolazione\mean_templates\Left-AG\4mean_baseline.ep');

confstruct.substract_column_at_start = 1; % detrend si o no
confstruct.method_GFPeak = 'GFPL2';       % gfp calculation method
confstruct.use_gfp_peaks = 1;             % utililla i picchi del gfp o tutto il file
confstruct.setall1 = 0;
confstruct.normalize = 0;
confstruct.tf = 512;
confstruct.fs = 256;
confstruct.minclusters = 1;
confstruct.maxclusters = 12;
confstruct.pmode = 0;
confstruct.algorithm = 1;
confstruct.similarity_measure = 'correlation';
confstruct.debug = 0;
confstruct.segment_min_len = 5;
confstruct.nch =27;

for k=1:size(elencoFiles,1)

    filename=elencoFiles(k).name;
    eeg = load(strcat(pathLeftAG_b,'\', filename));
    [ output ] = compute_mstate_parameters( confstruct, eeg, maps);
end
