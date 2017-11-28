function write_edf_Pierpaolo(file_out,DATA,header)

file_in = 'alonzi10_0101.edf'

%Open file and read header
ifp=fopen(file_in,'r');
ofp=fopen(file_out,'w');

vers=fscanf(ifp,'%c',8);  fprintf(ofp,'%c',vers);   % versione del formato
id=fscanf(ifp,'%c',80);   fprintf(ofp,'%c',id);      % identificazione del pz
rid=fscanf(ifp,'%c',80);  fprintf(ofp,'%c',rid);    % identificazione della registrazione
date=fscanf(ifp,'%c',8);  fprintf(ofp,'%c',date);
time=fscanf(ifp,'%c',8);  fprintf(ofp,'%c',time);
s_head=fscanf(ifp,'%c',8); fprintf(ofp,'%c',s_head);
res1=fscanf(ifp,'%c',44);  fprintf(ofp,'%c',res1);
nrec1=fscanf(ifp,'%c',8); 
reclen=fscanf(ifp,'%c',8); 

nrec_n = floor(length(DATA)/(str2num(reclen)*256));
nrec = '        ';
nrec(1:length(num2str(nrec_n))) = num2str(nrec_n);
 fprintf(ofp,'%c',nrec);
 fprintf(ofp,'%c',reclen);

nsig=fscanf(ifp,'%c',4);  fprintf(ofp,'%c',nsig);
nsig = sscanf(nsig,'%f');

for k=1:nsig
    if k<20
    LL(k,:)=fscanf(ifp,'%c',16);
    xx = '                ';
    xx(1:size(header.ch(k).label,2)) = header.ch(k).label;
     
    fprintf(ofp,'%c',xx);
    else
        LL(k,:)=fscanf(ifp,'%c',16);
    fprintf(ofp,'%c',LL(k,:));
    end
end
for k=1:nsig
    trasd=fscanf(ifp,'%c',80);
    fprintf(ofp,'%c',trasd);
end

for k=1:nsig
    phdim=fscanf(ifp,'%c',8);
    fprintf(ofp,'%c',phdim');
end
for k=1:nsig
    phmin=fscanf(ifp,'%c',8);
    fprintf(ofp,'%c',phmin);
end
for k=1:nsig
    phmax=fscanf(ifp,'%c',8);
    phmax_n(k) = sscanf(phmax,'%f');
    fprintf(ofp,'%c',phmax);
end
for k=1:nsig
    digmin=fscanf(ifp,'%c',8);
    fprintf(ofp,'%c',digmin);
end
for k=1:nsig
    digmax=fscanf(ifp,'%c',8);
    fprintf(ofp,'%c',digmax);
    digmax_n(k) = sscanf(digmax,'%f');
end

for k=1:nsig
    preft=fscanf(ifp,'%c',80);
    fprintf(ofp,'%c',preft');
end
for k=1:nsig
    nchpt=fscanf(ifp,'%c',8);
    fprintf(ofp,'%c',nchpt);
    Nchpt(k) = sscanf(nchpt,'%f');     %%%prova
end
for k=1:nsig
    reserved=fscanf(ifp,'%c',32);
    fprintf(ofp,'%c',reserved);
end

cal=phmax_n(1)./digmax_n(1);

n_ann = 0;
for k=1:nsig
    if strcmp(LL(k,:),'EDF Annotations ')
        n_ann = n_ann +1;
    end
end

for k=1:nsig-n_ann
    DATA(k,:) = DATA(k,:)/cal;
end

nrec = sscanf(nrec,'%f');
for sample = 1:nrec_n
    data_buffer = DATA(:,(sample-1)*Nchpt(1)+1:sample*Nchpt(1));
    for j=1:nsig-n_ann
        x = data_buffer(j,:);
        fwrite(ofp,x,'int16');
    end
    for j=n_ann:-1:1
        x = ones(1,Nchpt(end-j+1));
        fwrite(ofp,x,'int16');
    end
     
end
fclose('all')