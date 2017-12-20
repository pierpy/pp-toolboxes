function [elec, newFs, newtpdata] = itab_readedf(file, h_filtering, l_filtering, freq_notch, dec, tch)

    %Open file and read header
    disp('>> reading data...');

    ifp = fopen(file,'r');

    vers = fscanf(ifp,'%c',8);
    id = sscanf(fscanf(ifp,'%c',80),'%d');
    rid = fscanf(ifp,'%c',80);
    date = fscanf(ifp,'%c',8);
    time = fscanf(ifp,'%c',8);
    s_head = sscanf(fscanf(ifp,'%c',8),'%d');
    res1 = fscanf(ifp,'%c',44);
    nrec = sscanf(fscanf(ifp,'%c',8),'%f');
    reclen = fscanf(ifp,'%c',8);
    nsig = sscanf(fscanf(ifp,'%c',4),'%d');

    label = fscanf(ifp,'%c',[16,nsig])';
    trasd = fscanf(ifp,'%c',[80,nsig])';
    phdim = fscanf(ifp,'%c',[8,nsig])';
    phmin = sscanf(fscanf(ifp,'%c',[8,nsig]),'%f');
    phmax = sscanf(fscanf(ifp,'%c',[8,nsig]),'%f');
    digmin = sscanf(fscanf(ifp,'%c',[8,nsig]),'%f');
    digmax = sscanf(fscanf(ifp,'%c',[8,nsig]),'%f');
    preft = fscanf(ifp,'%c',[80,nsig])';
    nchpt = sscanf(fscanf(ifp,'%c',[8,nsig]),'%d');
    reserved = fscanf(ifp,'%c',[32,nsig])';

    c=1;
    cal = phmax./digmax;
    Fs = nchpt(1)/str2num(reclen);

    n_ann = 0;
    for k=1:nsig
        if strcmp(label(k,:),'EDF Annotations ')
            n_ann = n_ann +1;
        end
    end

    data = zeros(length(tch),nchpt(1)*nrec);

    for(sample=1:nrec)
        sample
        %Read data
        fseek(ifp,s_head+2*sum(nchpt)*(sample-1),'bof');

        meg = zeros(nchpt(1),nsig-n_ann);
        meg = fread(ifp,size(meg),'int16')';
        for(k = 1:nsig-n_ann)
            meg(k,:) = meg(k,:)*cal(k);
        end
        data(:,nchpt(1)*(sample-1)+1:nchpt(1)*(sample-1)+nchpt(1)) = meg(tch, :);
    end

    elec = zeros(length(tch),round(length(data)/dec));
    for cont = 1: length(tch);
        sig = data(cont,:);

        for z = 1: length(freq_notch)
            sig = itab_filter_notch(sig,freq_notch(z),Fs);
        end
        switch h_filtering
            case 0.5
                sig = itab_filter_h05(sig,Fs);
            case 1
                sig = itab_filter_h1(sig,Fs);
            case 40
                sig = itab_filter_h40(sig,Fs);
        end
        switch l_filtering
            case 40
                sig = itab_filter_l40(sig,Fs);
            case 80
                sig = itab_filter_l80(sig,Fs);
            case 100
                sig = itab_filter_l100(sig,Fs);
            case 150
                sig = itab_filter_l150(sig,Fs);
        end
        disp([ num2str(cont) ' / ' num2str(length(tch))]);

        sig1 = decimate(sig,dec);
        elec(cont,:) = sig1;
        clear sig sig1
    end
    disp(strcat('>> Loaded',num2str(cont),'channels'))
    disp(strcat('>> file name:', file(end-16:end-4)));
    disp(strcat('>> high pass frq: ', num2str(h_filtering), 'Hz'));
    disp(strcat('>> low pass frq: ', num2str(l_filtering), 'Hz'));
    disp(strcat('>> Notch frq: ', num2str(freq_notch), 'Hz'));
    disp(strcat('>> Decimator factor: ', num2str(dec)));
    newFs = round(Fs/dec);
    newtpdata = length(elec);
    disp(strcat('>> New sample frequency: ', num2str(newFs)));
    disp(strcat('>> New data leng: ', num2str(newtpdata)));
end

