tt = ones(1,230400);

% SKIPPED INTERVALS: ica_par.intervals.skipped_intervals
for k=1:size(ica_par.intervals.skipped_intervals,1)
    tt(ica_par.intervals.skipped_intervals(k,1):ica_par.intervals.skipped_intervals(k,2)) = 0;
end

best_iter = ICA.best_iter.index;


g_ic = ICA.iteration(best_iter).brain_ic;
A = ICA.iteration(best_iter).mixing;

data = A(:,g_ic)*ICA.best_iter.IC(g_ic,:);
dataArt = A*ICA.best_iter.IC;

% Le label IN: header.ch.label

ind_no = find(tt==0);


DATA = zeros(19,length(data)+length(ind_no));
npt = size(DATA,2);

tt(npt+1:end)=[];
DATA(:,find(tt==1))=data;

DATA(20,:)=0;

file_edf = 'pippetto.edf';
write_edf_Pierpaolo(file_edf,DATA,header)
