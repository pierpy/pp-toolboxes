it = ICA.best_iter.index;


good_ic = ICA.iteration(it). brain_ic;
XX = ICA.iteration(it).mixing(:,good_ic)*ICA.best_iter.IC(good_ic,:);
YY = ICA.iteration(it).mixing*ICA.best_iter.IC;
close all
figure
for k=1:19
    plot(1:length(YY),YY(k,:),'b',1:length(XX),XX(k,:),'r')
    pause
end

NN = 0;
tt = 1:1000000;
for k=1:size(ica_par.intervals.skipped_intervals,1)
    tt(ica_par.intervals.skipped_intervals(k,1):ica_par.intervals.skipped_intervals(k,2))=0;
    NN = NN + length(ica_par.intervals.skipped_intervals(k,1):ica_par.intervals.skipped_intervals(k,2));
end
L = length(XX)+NN;
tt = tt(1:L);
indX = find(tt==0);
tt(indX)=[];

if size(tt,2) ~= size(XX,2)
    error('NOOO')
end

close all
figure
plot(tt,XX)
cd('G:\DATI_Filippo\Stroke_acuti\MICROSTATI\FileTXT\Controlli')
nch = 19;
npt = header.smpfq*2;


ind =1;

[a1,b1] = butter(2,[2/(header.smpfq/2) 20/(header.smpfq/2)]);
[a2,b2] = butter(2,[1/(header.smpfq/2) 40/(header.smpfq/2)]);

for k=1:19
    XFilt1(k,:) = filtfilt(a1,b1,XX(k,:));
    XFilt2(k,:) = filtfilt(a2,b2,XX(k,:));
end


filename1 = 'valpie_F1_00.txt';
filename2 = 'valpie_F2_00.txt';
for k=1:30
    if k<10
        filename1(12)=num2str(k);
        filename2(12)=num2str(k);
    else
        filename1(11:12)=num2str(k);
        filename2(11:12)=num2str(k);
    end
    
    winS = ind+1:ind+npt;
    while ~isempty(find(diff(tt(winS))~=1))
        disp(strcat('sto qua',num2str(ind)))
        pause
        ind = ind+npt;
        winS = ind+1:ind+npt;
    end
    xx = XFilt1(:,winS);
    ifp = fopen(filename1,'w')
    for ii = 1:npt
        for jj = 1:nch
            fprintf(ifp,'%.4f\t',xx(jj,ii));
        end
        fprintf(ifp,'\n')
    end
    fclose all
    
     xx = XFilt2(:,winS);
    ifp = fopen(filename2,'w')
    for ii = 1:npt
        for jj = 1:nch
            fprintf(ifp,'%.4f\t',xx(jj,ii));
        end
        fprintf(ifp,'\n')
    end
    fclose all
    
    ind = ind+npt;
end

    
    