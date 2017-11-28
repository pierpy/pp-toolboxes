function RaguDumpToVision(rd)

DirName = uigetdir('','Select directory to store the data');
if DirName == 0
    return
end

for s = 1:size(rd.Names,1)
    for c = 1:size(rd.Names,2)
        fn = rd.Names{s,c};
        idx = find(fn == '.',1,'last');
        if isempty(idx)
            idx = numel(fn)+1;
        end
        OutFile = fn(1:(idx-1));
        out = squeeze(rd.V(s,c,:,:));
        outt = out';
        save([DirName '\\' OutFile '.asc'],'outt','-ASCII');

        [fid,Message] = fopen([DirName '\\' OutFile '.vhdr'],'wt');

        if fid == -1
            error(Message);
        end

        fprintf(fid,'Brain Vision Data Exchange Header File Version 1.0\n');
        fprintf(fid,'[Common Infos]\n');
        fprintf(fid,['DataFile=' OutFile '.asc\n']);
        fprintf(fid,'DataFormat=ASCII\n');
        fprintf(fid,'DataOrientation=MULTIPLEXED\n');
        fprintf(fid,'DataType=TIMEDOMAIN\n');
        fprintf(fid,'NumberOfChannels=%i\n',size(out,1));
        fprintf(fid,'DataPoints=%i\n',size(out,2));
        fprintf(fid,'Averaged=YES\n');
        fprintf(fid,'; Sampling interval in microseconds if time domain (convert to Hertz:\n');
        fprintf(fid,'; 1000000 / SamplingInterval) or in Hertz if frequency domain:\n');
        fprintf(fid,'SamplingInterval=%i\n',1000*rd.DeltaX);
        fprintf(fid,'[Channel Infos]\n');
 
        for i = 1:size(out,1)
            fprintf(fid,'Ch%i=%s,,1\n',i,rd.Channel(i).Name);
        end

        fprintf(fid,'\n[Coordinates]\n');
        
        if isfield(rd.Channel(1),'Radius')
            for i = 1:size(out,1)
                fprintf(fid,'Ch%i=%f,%f,%f\n',i,rd.Channel(i).Radius,rd.Channel(i).Theta,rd.Channel(i).Phi);
            end
        else   
            for i = 1:size(out,1)
                fprintf(fid,'Ch%i=%f,%f,%f\n',i,rd.Channel(i).CoordRadius,rd.Channel(i).CoordTheta,rd.Channel(i).CoordPhi);
            end
        end
        fclose(fid);
        
        
        [fid,Message] = fopen([DirName '\\' OutFile '.vmrk'],'wt');

        if fid == -1
            error(Message);
        end

        fprintf(fid,'Brain Vision Data Exchange Marker File, Version 1.0\n');
        fprintf(fid,'[Common Infos]\n');
        fprintf(fid,['DataFile=' OutFile '.asc\n']);

        fprintf(fid,'[Marker Infos]\n');
        fprintf(fid,'Mk1=New Segment,,1,1,0,00000000000000000000\n');
        fprintf(fid,'Mk2=Time 0,,%i,1,0\n',round(-rd.TimeOnset * rd.DeltaX)+1);
        fclose(fid);
    end
end
     


