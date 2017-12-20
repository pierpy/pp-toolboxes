function [fullpathparfile] = itab_writeparameters(subjfolderpath, filename)
    parfilename = strcat(filename, '.par');
    fullpathparfile = fullfile(subjfolderpath, parfilename);
    fp = fopen(fullpathparfile,'wt'); 

    risposta = inputdlg({'High-pass filtering (0.5, 1 or 40 Hz)','Low-pass filtering (40, 80, 100 or 150 Hz)','Frequencies to be excluded (in Hz)',...
        'Decimation factor','FastICA nonlinearity (tanh, gauss, both)',...
        'Number of segments to be excluded'},'Processing parameter selection',1,{'1','150','none','3','tanh','0','0'});
    h_filtering = str2num(risposta{1});
    l_filtering = str2num(risposta{2});
    freq_notch = str2num(risposta{3});
    dec=str2num(risposta{4});
    ica_method = risposta{5};
    num_segments = str2num(risposta{6});

    if not( h_filtering==0.5 | h_filtering==1 | h_filtering==40)
        h_filtering=1;
    end
    if not(l_filtering == 40 | l_filtering == 80 | l_filtering == 100 | l_filtering == 150)
        l_filtering=150;
    end

    ica_iter = 5;

    fprintf(fp, '%s\n\n', 'ICA PARAMETER SETTINGS');

    fprintf(fp, '%s %s\n', 'High-pass filter cut off frequency [Hz]: ', num2str(h_filtering));
    fprintf(fp, '%s %s\n', 'Low-pass filter cut off frequency [Hz]: ', num2str(l_filtering));
    fprintf(fp, '%s %s\n', 'Notch filter cut off frequencies [Hz]: ', risposta{3});
    fprintf(fp, '%s %s\n', 'Decimation factor: ', risposta{4});
    fprintf(fp, '%s %s\n', 'ICA method: ', risposta{5});
    fprintf(fp, '%s %s\n', 'ICA iteration number: ',num2str(ica_iter));
    fprintf(fp, '%s %s\n', 'Number of segments to be excluded:', risposta{6});

    if(num_segments~=0)
        fprintf(fp, '%s \n', 'Skipped Intervals [s]:');

        for z=1:num_segments
            risposta_new = inputdlg({'Interval onset [s]','Interval offset [s]'},['Data segment to be excluded - no.' num2str(z)],1,{'',''});
            fprintf(fp, '\t%.2f \t%.2f\n', str2num(risposta_new{1}), str2num(risposta_new{2}));
        end
    end

    fprintf(fp, '\n%s \t%s\n', 'Study type:', 'ongoing');

    tch=[125,126,127,128];
    nfigs=length(tch);
    risposta = inputdlg({'Fill in the number of electric channels'},'Electric channel selection for artifact detection',1,{num2str([0:nfigs-1])});
    ech = [str2num(risposta{1})+1];
    close all

    fprintf(fp, '\n%s\n','Selected electric channels');
    fprintf(fp, '%d\t',tch(ech));  
    close;

    ButtonName=questdlg('Do you want to load respiration channel?', ...
    'Question', ...
        'Yes','No','No');
    fprintf(fp, '\n%s\t%s\n','Respiration channel loaded:', ButtonName);  

    fclose(fp);
end
